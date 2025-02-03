import torch,sys,time,tqdm
from transformers import AutoModel, AutoModelForMaskedLM, AutoTokenizer
from torch.utils.data import Dataset, DataLoader
from Bio import SeqIO
import numpy as np
import pandas as pd
import pickle
import argparse, sys, os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", dest="inputFASTA", type=str, default=None,
                        help="The directory of input fasta")
    parser.add_argument("-output", dest = "output", help = "The directory of output")
    parser.add_argument("-model", dest = "modelDir", default="facebook/esm1b_t33_650M_UR50S", help = "The directory of pre-trained model")
    parser.add_argument("-device", dest = "device", default="cuda:0", help = "The device to run the model")
    args = parser.parse_args()
    return args

def get_model(model_path, device = 'cpu'):
    tokenizer = AutoTokenizer.from_pretrained(model_path)
    model = AutoModelForMaskedLM.from_pretrained(model_path)
    model.to(device)
    model.eval()
    return model, tokenizer


class SequenceDataset(Dataset):
    def __init__(self, sequences, tokenizer, names):
        self.sequences = sequences
        self.tokenizer = tokenizer
        self.names = names

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        name = self.names[idx]
        original_length = len(sequence)
        chunks = []

        if original_length > 1022:
            for i in range(0, original_length, 511):
                chunk = sequence[i:i + 1022]
                encoded_sequence = self.tokenizer.encode_plus(
                    chunk,
                    add_special_tokens=True,
                    truncation=True,
                    max_length=1024,
                    padding="max_length",
                    return_attention_mask=True,
                    return_tensors="pt"
                )
                chunks.append({
                    'sequence': chunk,
                    'input_ids': encoded_sequence['input_ids'],
                    'attention_mask': encoded_sequence['attention_mask'],
                    'seqlen': len(chunk),
                    'name': name
                })
                if i + 1022 >= original_length:
                    break
        else:
            encoded_sequence = self.tokenizer.encode_plus(
                sequence,
                add_special_tokens=True,
                truncation=True,
                max_length=original_length + 2,
                padding="max_length",
                return_attention_mask=True,
                return_tensors="pt"
            )
            chunks.append({
                'sequence': sequence,
                'input_ids': encoded_sequence['input_ids'],
                'attention_mask': encoded_sequence['attention_mask'],
                'seqlen': original_length,
                'name': name
            })

        return chunks


args = parse_args()
device = args.device
model, tokenizer = get_model(args.modelDir, device = device)
sequences = []
names = []
seqLen = []
for record in SeqIO.parse(args.inputFASTA, "fasta"):
    sequences.append(str(record.seq))
    names.append(record.id)
    seqLen.append(len(str(record.seq)))

dataset = SequenceDataset(sequences = sequences, names = names,
                          tokenizer=tokenizer)
loader = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=1)
resDic = {}
amino_acids = list('LAGVSERTIDPKQNFYMHWC')
for batch in tqdm.tqdm(loader, desc="ESM logits...", unit="batch"):
    if len(batch) > 1:
        probsRes = []
        for idx, chunk in enumerate(batch):
            input_ids = chunk['input_ids'].to(device).squeeze(1)
            batch_len = chunk['seqlen'].numpy()
            transcriptID = chunk['name']
            attention_mask = chunk['attention_mask'].to(device).squeeze(1)
            with torch.inference_mode():
                all_logits = model(input_ids=input_ids,attention_mask=attention_mask).logits
            target_logits = all_logits[0, 1:(batch_len[0]+1), [tokenizer.get_vocab()[aa] for aa in amino_acids]]
            probs = torch.nn.functional.softmax(target_logits.cpu(), dim=1).numpy()
            if idx == 0:
                probs = probs
            else:
                probs = probs[511:]

            probsRes.append(probs)
            probs = np.vstack(probsRes)
    else:
        batch = batch[0]
        input_ids = batch['input_ids'].to(device).squeeze(1)
        batch_len = batch['seqlen'].numpy()
        transcriptID = batch['name']
        attention_mask = batch['attention_mask'].to(device).squeeze(1)
        with torch.inference_mode():
            all_logits = model(input_ids=input_ids,attention_mask=attention_mask).logits
        target_logits = all_logits[0, 1:(batch_len[0]+1), [tokenizer.get_vocab()[aa] for aa in amino_acids]]
        probs = torch.nn.functional.softmax(target_logits.cpu(), dim=1).numpy()
    resDic[transcriptID[0]] = probs

output = args.output

np.savez_compressed(output, **resDic)