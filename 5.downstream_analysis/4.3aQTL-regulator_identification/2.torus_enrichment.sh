#!/bin/bash
ct=$1
torus -d  01.${ct}.assoc.txt.gz  -smap 01.${ct}.snpmap.txt.gz -gmap 01.${ct}.genemap.txt.gz -annot est/01.Exc.snpanno.txt.gz  -est > est_res/${ct}.snpanno.est

