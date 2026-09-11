import torch,sys,time,tqdm
from transformers import AutoModelForMaskedLM, AutoTokenizer
from torch.utils.data import Dataset, DataLoader
from Bio import SeqIO
import numpy as np
import pandas as pd
import pickle
import argparse, sys, os

# PlantCaduceus (kuleshov-group/PlantCaduceus_l32) is a DNA-sequence Caduceus/Mamba model, not
# a standard transformer - same overall extraction structure as src/4_ESM_logits.py, but with
# a 4-letter lowercase nucleotide vocab instead of the 20-amino-acid one, no leading/trailing
# special token added by the tokenizer (confirmed empirically: encoding an N-base sequence
# yields exactly N token ids, unlike ESM's leading <cls>), and a custom model class requiring
# `trust_remote_code=True`. Requires a CUDA GPU - PlantCaduceus's mamba_ssm backend has no CPU
# fallback (confirmed: a CPU forward pass fails inside its Triton fused-layernorm kernel with
# "invalid argument to exchangeDevice"), matching this stage's existing SCINET-GPU-only
# convention for LLM scoring.

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", dest="inputFASTA", type=str, default=None,
                        help="The directory of input fasta (nucleotide sequences)")
    parser.add_argument("-output", dest = "output", help = "The directory of output")
    parser.add_argument("-model", dest = "modelDir", default="kuleshov-group/PlantCaduceus_l32", help = "The directory of pre-trained model")
    parser.add_argument("-device", dest = "device", default="cuda:0", help = "The device to run the model")
    parser.add_argument("-max_length", dest = "maxLength", type=int, default=4096,
                        help = "Chunk sequences longer than this many bases (overlapping, matching 4_ESM_logits.py's approach)")
    args = parser.parse_args()
    return args

def get_model(model_path, device = 'cpu'):
    tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(model_path, trust_remote_code=True)
    model.to(device)
    model.eval()
    return model, tokenizer


class SequenceDataset(Dataset):
    def __init__(self, sequences, tokenizer, names, max_length):
        self.sequences = sequences
        self.tokenizer = tokenizer
        self.names = names
        self.max_length = max_length

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        name = self.names[idx]
        original_length = len(sequence)
        half_chunk = self.max_length // 2
        chunks = []

        if original_length > self.max_length:
            for i in range(0, original_length, half_chunk):
                chunk = sequence[i:i + self.max_length]
                encoded_sequence = self.tokenizer.encode_plus(
                    chunk,
                    add_special_tokens=True,
                    truncation=True,
                    max_length=self.max_length,
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
                if i + self.max_length >= original_length:
                    break
        else:
            encoded_sequence = self.tokenizer.encode_plus(
                sequence,
                add_special_tokens=True,
                truncation=True,
                max_length=original_length,
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
    sequences.append(str(record.seq).lower())
    names.append(record.id)
    seqLen.append(len(str(record.seq)))

dataset = SequenceDataset(sequences = sequences, names = names,
                          tokenizer=tokenizer, max_length=args.maxLength)
loader = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=1)
resDic = {}
nucleotides = list('acgt')
half_chunk = args.maxLength // 2
for batch in tqdm.tqdm(loader, desc="PlantCAD logits...", unit="batch"):
    if len(batch) > 1:
        probsRes = []
        for idx, chunk in enumerate(batch):
            input_ids = chunk['input_ids'].to(device).squeeze(1)
            batch_len = chunk['seqlen'].numpy()
            transcriptID = chunk['name']
            attention_mask = chunk['attention_mask'].to(device).squeeze(1)
            with torch.inference_mode():
                all_logits = model(input_ids=input_ids,attention_mask=attention_mask).logits
            # no leading special token (unlike ESM's <cls>) - positions are 0-indexed
            target_logits = all_logits[0, 0:batch_len[0], [tokenizer.get_vocab()[nt] for nt in nucleotides]]
            probs = torch.nn.functional.softmax(target_logits.cpu(), dim=1).numpy()
            if idx == 0:
                probs = probs
            else:
                probs = probs[half_chunk:]

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
        target_logits = all_logits[0, 0:batch_len[0], [tokenizer.get_vocab()[nt] for nt in nucleotides]]
        probs = torch.nn.functional.softmax(target_logits.cpu(), dim=1).numpy()
    resDic[transcriptID[0]] = probs

output = args.output

np.savez_compressed(output, **resDic)
