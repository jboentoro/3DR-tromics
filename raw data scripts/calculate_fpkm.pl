opendir(DIR1,"2_map");
while($files =readdir(DIR1)){
if($files =~ /\.summ$/){

open(IN1,"2_map/$files")|| die "cant find 1\n";
while(<IN1>){
chomp $_;
if($_ =~ /aligned concordantly exactly 1 time/){
@arr1=split('\(',$_);
$arr1[0]=~ s/\s+//g;  #unique count
}
}

$file2=$files;
$file3=$files;
$file4=$files;
$file2 =~ s/\.summ$/\.unq\.cov/g;
$file3 =~ s/\.summ$/\.fpkm\.txt/g;
$file4 =~ s/\.summ$//g;

open(OUT1,">2_map/$file3");
print OUT1 "geneid\t",$file4,"\n";
open(IN2,"2_map/$file2");
while(<IN2>){
chomp $_;
@arr2=split('\t',$_);
@arr3=split(';',$arr2[8]);
$arr3[0] =~ s/ID\=//g;
if($arr2[10] != 0){$fpkmcal=(($arr2[9]/$arr2[10])/$arr1[0])*1000000000;}
else{$fpkmcal="NA";}
print OUT1 $arr3[0],"\t",$fpkmcal,"\n";
}
close IN1;
close OUT1;
close IN2;
}
}

