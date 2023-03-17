perl -ane '
    open $out{$F[0]}, ">", $F[0]."_junctions_1.bed" unless $out{$F[0]}; 
    print { $out{$F[0]} } $_;
' junctions_1.bed
