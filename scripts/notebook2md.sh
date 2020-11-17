#!/bin/bash
for notebook in notebooks/*.ipynb; do
   jupyter nbconvert\
      --to markdown\
      --output-dir docs\
      --TemplateExporter.exclude_input=True\
      $notebook
done

