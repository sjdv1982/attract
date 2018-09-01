#!/bin/bash
d=`dirname "$0"`
for i in `awk '{print $1}' $d/list_lib_gag`; do wget http://www.prot-gag.org/uploads/$i.lib -O $d/$i.lib; done
