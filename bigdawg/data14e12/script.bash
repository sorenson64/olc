#!/bin/bash
for i in *.out
do
  tail -n 1 "$i" >> bench.txt
done

