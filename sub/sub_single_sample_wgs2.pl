#!/usr/bin/perl -w
use strict;

###USER INPUT PARAMETERS###
my $name="a.single_sample2";
my $par=2;
my $script="single_sample_wgs_pipeline2.pl";
my $opt="sc";
my $i="318089";
my $j="SEQ 1 40";
my $k="NA";
my $node=1;
my $ppn=1;
my $mem="20g";
my $scratch="20g";
my $dir="/hpf/projects/sbowdin/andyding/yiming";
#my $wkt="359:59:59";
my $wkt="23:59:59";
#my $wkt="1:59:59";
my @module=("R");
###########################

# version 20220331 - update this date if anything below changes
my @i;
my @split1=split /\s/,$i;
if ($split1[0] eq "SEQ"){
 for (my $i=$split1[1]; $i<=$split1[2]; $i++){
  push @i, $i;
 }
}
elsif ($split1[0] eq "MORE"){
 open file1, "<$split1[1]";
 while (<file1>){
  chomp;
  my @split2=split /\s/,$_;
  foreach my $split2 (@split2){
   push @i, $split2;
  }
 }
}
else{
 foreach my $split1 (@split1){
  push @i, $split1;
 }
}

my @j;
@split1=split /\s/,$j;
if ($split1[0] eq "SEQ"){
 for (my $i=$split1[1]; $i<=$split1[2]; $i++){
  push @j, $i;
 }
}
elsif ($split1[0] eq "MORE"){
 open file1, "<$split1[1]";
 while (<file1>){
  chomp;
  my @split2=split /\s/,$_;
  foreach my $split2 (@split2){
   push @j, $split2;
  }
 }
}
else{
 foreach my $split1 (@split1){
  push @j, $split1;
 }
}

my @k;
@split1=split /\s/,$k;
if ($split1[0] eq "SEQ"){
 for (my $i=$split1[1]; $i<=$split1[2]; $i++){
  push @k, $i;
 }
}
elsif ($split1[0] eq "MORE"){
 open file1, "<$split1[1]";
 while (<file1>){
  chomp;
  my @split2=split /\s/,$_;
  foreach my $split2 (@split2){
   push @k, $split2;
  }
 }
}
else{
 foreach my $split1 (@split1){
  push @k, $split1;
 }
}

if ($par==1){
 foreach my $i (@i){
  open out1, ">./$name.$i";
  print out1 "#!/bin/bash\n";
  print out1 "#PBS -l walltime=$wkt\n";
  print out1 "#PBS -l nodes=$node:ppn=$ppn\n";
  print out1 "#PBS -l mem=$mem,vmem=$mem,file=$scratch\n";
  foreach my $module (@module){
   print out1 "module load $module\n";
  }
  print out1 "cd $dir\n";
  print out1 "perl $dir/$script -$opt $i\n";
  close out1;
  system ("qsub $name.$i");
  system ("rm $name.$i");
 }
}

if ($par==2){
 my @split1=split //,$opt;
 foreach my $i (@i){
  foreach my $j (@j){
   open out1, ">./$name.$i.$j";
   print out1 "#!/bin/bash\n";
   print out1 "#PBS -l walltime=$wkt\n";
   print out1 "#PBS -l nodes=$node:ppn=$ppn\n";
   print out1 "#PBS -l mem=$mem,vmem=$mem,file=$scratch\n";
   foreach my $module (@module){
    print out1 "module load $module\n";
   }
   print out1 "cd $dir\n";
   print out1 "perl $dir/$script -$split1[0] $i -$split1[1] $j\n";
   close out1;
   system ("qsub $name.$i.$j");
   system ("rm $name.$i.$j");
  }
 }
}

if ($par==3){
 my @split1=split //,$opt;
 foreach my $i (@i){
  foreach my $j (@j){
   foreach my $k (@k){
    open out1, ">./$name.$i.$j.$k";
    print out1 "#!/bin/bash\n";
    print out1 "#PBS -l walltime=$wkt\n";
    print out1 "#PBS -l nodes=$node:ppn=$ppn\n";
    print out1 "#PBS -l mem=$mem,vmem=$mem,file=$scratch\n";
    foreach my $module (@module){
     print out1 "module load $module\n";
    }
    print out1 "cd $dir\n";
    print out1 "perl $dir/$script -$split1[0] $i -$split1[1] $j -$split1[2] $k\n";
    close out1;
    system ("qsub $name.$i.$j.$k");
    system ("rm $name.$i.$j.$k");
   }
  }
 }
}

exit 2;
