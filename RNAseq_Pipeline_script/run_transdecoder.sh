#!/usr/bin/env bash

transcripts=$1

~/src/trinityrnaseq_r2013-02-25/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl -t $transcripts -S
