@rem
@rem Copyright 2015 the original author or authors.
@rem
@rem Licensed under the Apache License, Version 2.0 (the "License");
@rem you may not use this file except in compliance with the License.
@rem You may obtain a copy of the License at
@rem
@rem      https://www.apache.org/licenses/LICENSE-2.0
@rem
@rem Unless required by applicable law or agreed to in writing, software
@rem distributed under the License is distributed on an "AS IS" BASIS,
@rem WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
@rem See the License for the specific language governing permissions and
@rem limitations under the License.
@rem

@if "%DEBUG%"=="" @echo off
@rem ##########################################################################
@rem
@rem  patch-scaffolds startup script for Windows
@rem
@rem ##########################################################################

@rem Set local scope for the variables with windows NT shell
if "%OS%"=="Windows_NT" setlocal

set DIRNAME=%~dp0
if "%DIRNAME%"=="" set DIRNAME=.
@rem This is normally unused
set APP_BASE_NAME=%~n0
set APP_HOME=%DIRNAME%..

@rem Resolve any "." and ".." in APP_HOME to make it shorter.
for %%i in ("%APP_HOME%") do set APP_HOME=%%~fi

@rem Add default JVM options here. You can also use JAVA_OPTS and PATCH_SCAFFOLDS_OPTS to pass JVM options to this script.
set DEFAULT_JVM_OPTS=

@rem Find java.exe
if defined JAVA_HOME goto findJavaFromJavaHome

set JAVA_EXE=java.exe
%JAVA_EXE% -version >NUL 2>&1
if %ERRORLEVEL% equ 0 goto execute

echo.
echo ERROR: JAVA_HOME is not set and no 'java' command could be found in your PATH.
echo.
echo Please set the JAVA_HOME variable in your environment to match the
echo location of your Java installation.

goto fail

:findJavaFromJavaHome
set JAVA_HOME=%JAVA_HOME:"=%
set JAVA_EXE=%JAVA_HOME%/bin/java.exe

if exist "%JAVA_EXE%" goto execute

echo.
echo ERROR: JAVA_HOME is set to an invalid directory: %JAVA_HOME%
echo.
echo Please set the JAVA_HOME variable in your environment to match the
echo location of your Java installation.

goto fail

:execute
@rem Setup the command line

set CLASSPATH=%APP_HOME%\lib\patch-scaffolds.jar;%APP_HOME%\lib\biokotlin-0.15.jar;%APP_HOME%\lib\krangl-0.18.jar;%APP_HOME%\lib\dataframe-0.8.0-rc-7.jar;%APP_HOME%\lib\klaxon-5.6.jar;%APP_HOME%\lib\kotlin-reflect-1.9.10.jar;%APP_HOME%\lib\kotlin-statistics-1.2.1.jar;%APP_HOME%\lib\kotlin-logging-jvm-5.0.0.jar;%APP_HOME%\lib\ktor-client-cio-jvm-1.6.3.jar;%APP_HOME%\lib\ktor-client-core-jvm-1.6.3.jar;%APP_HOME%\lib\ktor-http-cio-jvm-1.6.3.jar;%APP_HOME%\lib\ktor-http-jvm-1.6.3.jar;%APP_HOME%\lib\ktor-network-tls-jvm-1.6.3.jar;%APP_HOME%\lib\ktor-network-jvm-1.6.3.jar;%APP_HOME%\lib\ktor-utils-jvm-1.6.3.jar;%APP_HOME%\lib\ktor-io-jvm-1.6.3.jar;%APP_HOME%\lib\kotlinx-coroutines-core-jvm-1.5.2.jar;%APP_HOME%\lib\plot-api-jvm-3.2.0.jar;%APP_HOME%\lib\lets-plot-common-2.3.0.jar;%APP_HOME%\lib\clikt-jvm.jar;%APP_HOME%\lib\mordant-jvm.jar;%APP_HOME%\lib\plot-config-portable-jvm-2.3.0.jar;%APP_HOME%\lib\plot-builder-portable-jvm-2.3.0.jar;%APP_HOME%\lib\plot-base-portable-jvm-2.3.0.jar;%APP_HOME%\lib\plot-common-portable-jvm-2.3.0.jar;%APP_HOME%\lib\vis-svg-portable-jvm-2.3.0.jar;%APP_HOME%\lib\base-portable-jvm-2.3.0.jar;%APP_HOME%\lib\colormath-jvm.jar;%APP_HOME%\lib\markdown-jvm-0.5.2.jar;%APP_HOME%\lib\kotlin-stdlib-jdk8-1.9.10.jar;%APP_HOME%\lib\kotlinx-serialization-json-jvm-1.2.2.jar;%APP_HOME%\lib\fuel-2.3.1.jar;%APP_HOME%\lib\kotlinx-datetime-jvm-0.3.1.jar;%APP_HOME%\lib\kotlinx-serialization-core-jvm-1.2.2.jar;%APP_HOME%\lib\kotlin-logging-jvm-2.0.5.jar;%APP_HOME%\lib\kotlin-stdlib-jdk7-1.9.10.jar;%APP_HOME%\lib\result-3.1.0.jar;%APP_HOME%\lib\kotlin-stdlib-1.9.22.jar;%APP_HOME%\lib\log4j-core-2.23.1.jar;%APP_HOME%\lib\poi-ooxml-5.2.2.jar;%APP_HOME%\lib\poi-5.2.2.jar;%APP_HOME%\lib\poi-ooxml-lite-5.2.2.jar;%APP_HOME%\lib\xmlbeans-5.0.3.jar;%APP_HOME%\lib\log4j-api-2.23.1.jar;%APP_HOME%\lib\htsjdk-4.1.0.jar;%APP_HOME%\lib\annotations-13.0.jar;%APP_HOME%\lib\kotlin-script-runtime-1.9.10.jar;%APP_HOME%\lib\graal-sdk-21.2.0.jar;%APP_HOME%\lib\commons-csv-1.8.jar;%APP_HOME%\lib\guava-30.1.1-jre.jar;%APP_HOME%\lib\gremlin-core-3.5.1.jar;%APP_HOME%\lib\jgrapht-core-1.5.1.jar;%APP_HOME%\lib\logback-classic-1.2.6.jar;%APP_HOME%\lib\fastutil-8.5.12.jar;%APP_HOME%\lib\lz4-java-1.8.0.jar;%APP_HOME%\lib\commons-jexl-2.1.1.jar;%APP_HOME%\lib\commons-configuration2-2.7.jar;%APP_HOME%\lib\commons-beanutils-1.9.4.jar;%APP_HOME%\lib\commons-logging-1.2.jar;%APP_HOME%\lib\snappy-java-1.1.10.5.jar;%APP_HOME%\lib\commons-compress-1.24.0.jar;%APP_HOME%\lib\xz-1.9.jar;%APP_HOME%\lib\json-20230618.jar;%APP_HOME%\lib\nashorn-core-15.4.jar;%APP_HOME%\lib\ngs-java-2.9.0.jar;%APP_HOME%\lib\commons-math3-3.6.1.jar;%APP_HOME%\lib\kotlin-jupyter-api-annotations-0.10.0-131-1.jar;%APP_HOME%\lib\arrow-memory-netty-8.0.0.jar;%APP_HOME%\lib\arrow-vector-8.0.0.jar;%APP_HOME%\lib\failureaccess-1.0.1.jar;%APP_HOME%\lib\listenablefuture-9999.0-empty-to-avoid-conflict-with-guava.jar;%APP_HOME%\lib\arrow-memory-core-8.0.0.jar;%APP_HOME%\lib\jsr305-3.0.2.jar;%APP_HOME%\lib\checker-qual-3.8.0.jar;%APP_HOME%\lib\error_prone_annotations-2.5.1.jar;%APP_HOME%\lib\j2objc-annotations-1.3.jar;%APP_HOME%\lib\gremlin-shaded-3.5.1.jar;%APP_HOME%\lib\commons-collections-3.2.2.jar;%APP_HOME%\lib\commons-text-1.9.jar;%APP_HOME%\lib\commons-lang3-3.11.jar;%APP_HOME%\lib\snakeyaml-1.27.jar;%APP_HOME%\lib\javatuples-1.2.jar;%APP_HOME%\lib\hppc-0.7.1.jar;%APP_HOME%\lib\jcabi-manifests-1.1.jar;%APP_HOME%\lib\javapoet-1.8.0.jar;%APP_HOME%\lib\exp4j-0.4.8.jar;%APP_HOME%\lib\jcl-over-slf4j-1.7.25.jar;%APP_HOME%\lib\slf4j-api-1.7.32.jar;%APP_HOME%\lib\jheaps-0.13.jar;%APP_HOME%\lib\logback-core-1.2.6.jar;%APP_HOME%\lib\asm-commons-7.3.1.jar;%APP_HOME%\lib\asm-util-7.3.1.jar;%APP_HOME%\lib\asm-analysis-7.3.1.jar;%APP_HOME%\lib\asm-tree-7.3.1.jar;%APP_HOME%\lib\asm-7.3.1.jar;%APP_HOME%\lib\netty-buffer-4.1.72.Final.jar;%APP_HOME%\lib\netty-common-4.1.72.Final.jar;%APP_HOME%\lib\arrow-format-8.0.0.jar;%APP_HOME%\lib\jackson-annotations-2.13.2.jar;%APP_HOME%\lib\jackson-datatype-jsr310-2.13.2.jar;%APP_HOME%\lib\jackson-databind-2.13.2.2.jar;%APP_HOME%\lib\jackson-core-2.13.2.jar;%APP_HOME%\lib\commons-codec-1.15.jar;%APP_HOME%\lib\flatbuffers-java-1.12.0.jar;%APP_HOME%\lib\commons-io-2.11.0.jar;%APP_HOME%\lib\curvesapi-1.07.jar;%APP_HOME%\lib\commons-collections4-4.4.jar;%APP_HOME%\lib\jcabi-log-0.14.jar;%APP_HOME%\lib\SparseBitSet-1.2.jar;%APP_HOME%\lib\jna-5.13.0.jar;%APP_HOME%\lib\fastutil-core-8.5.12.jar


@rem Execute patch-scaffolds
"%JAVA_EXE%" %DEFAULT_JVM_OPTS% %JAVA_OPTS% %PATCH_SCAFFOLDS_OPTS%  -classpath "%CLASSPATH%" net.maizegenetics.patchScaffolds.cli.PatchScaffoldsKt %*

:end
@rem End local scope for the variables with windows NT shell
if %ERRORLEVEL% equ 0 goto mainEnd

:fail
rem Set variable PATCH_SCAFFOLDS_EXIT_CONSOLE if you need the _script_ return code instead of
rem the _cmd.exe /c_ return code!
set EXIT_CODE=%ERRORLEVEL%
if %EXIT_CODE% equ 0 set EXIT_CODE=1
if not ""=="%PATCH_SCAFFOLDS_EXIT_CONSOLE%" exit %EXIT_CODE%
exit /b %EXIT_CODE%

:mainEnd
if "%OS%"=="Windows_NT" endlocal

:omega
