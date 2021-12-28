BEGIN{ FS="\t";OFS="\t" }
{ if(!NF || /^$/) { next } }
{ match($1,/:([0-9]+)-/,offset); $4=$4+offset[1]-1; $5=$5+offset[1]-1; gsub(/_sliding:[0-9]+-[0-9]+/,"",$1); print $0 }