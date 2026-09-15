@echo off
REM Usage examples:
REM   run_GRNBoost2.bat normal mesophyll
REM   run_GRNBoost2.bat drought mesophyll

set CONDITION=%1
set CELLTYPE=%2

call conda activate scenic

cd data\%CONDITION%\

pyscenic grn ^
  %CELLTYPE%_%CONDITION%_expression.tsv ^
  %CELLTYPE%_%CONDITION%_TF_list.txt ^
  -o %CELLTYPE%_%CONDITION%_GRNBoost2.tsv ^
  --method grnboost2 ^
  --num_workers 4
