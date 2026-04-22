#!/usr/bin/perl
use Inline::Files;
use utf8;
use open ":encoding(gbk)",":std";

#不能读入FILE请安装cpanm install Inline::Files
#如果不能读入数据，请将文件改为不带中文的名称

$t=0;

sub t1
{
    my $num=$_[0]%2;#取余
    if ($num==0)
    {
        print("偶数\n");
    }
    else
    {
        print("奇数\n");
    }
}

sub t2
{
    my @dir=glob($_[0].'/*');
    foreach my $filename(@dir)
    {
        if (-f $filename)
        {
            print{OUT}("$filename\n");
        }
    }
    foreach my $filename(@dir)
    {
        if (-d $filename)
        {
            t2($filename);
        }
    }
}

sub t5_1
{
    foreach my $n(0..$#_)
    {
        $_[$n]=~tr/ATCG/TAGC/;
        print("$_[$n]\n");
    }
}

sub t5_2
{
    my $count=0;
    my $length=0;
    foreach my $n(0..$#_)
    {
        $length=$length+length($_[$n]);
        if ($_[$n]=~/CG/g)
        {
            $count++;
        }
        if ($_[$n]=~/GC/g)
        {
            $count++;
        }
    }
    #print("$count $length\n");
    print("比例为:");
    printf("%f%\n",$count/($length-1)*100);
}

sub t6
{
    my @length=();
    my $result=0;
    for(my $i=0;$i<$count-1;$i++)
    {
        for(my $j=$i+1;$j<$count-1;$j++)
        {
            $result=(($_[$i][2]-$_[$j][2])**2+($_[$i][3]-$_[$j][3])**2+($_[$i][4]-$_[$j][4])**2)**1/3;
            push(@length,[$_[$i][0],$_[$i][1],$_[$j][0],$_[$j][1],$result]);
        }
    }
    foreach (0..$#length)
    {
        print("@{$length[$_]}\n");
    }
    my @max=sort{$length[$b][4]<=>$length[$a][4]} (0..$#length);
    print("最大值为:");
    print("@{$length[$max[0]]}\n");
}

#1
$t++;
print("-"x20,"$t","-"x20,"\n");
chomp($NUM=<STDIN>);
t1($NUM);
#2
$t++;
print("-"x20,"$t","-"x20,"\n");
open OUT,">","filepath.txt";
print("请输入文件路径,如C:(尾部省略'/'):");
chomp($path=<STDIN>);
t2($path);
close OUT;
#3
$t++;
print("-"x20,"$t","-"x20,"\n");
chomp($array=<STDIN>);
if ($array=~/[0-9a-z]+@[0-9a-z]+\.com/i)
{
    print("合法\n");
}
else{
    print("不合法\n");
}
#4
$t++;
%hash=();
print("-"x20,"$t","-"x20,"\n");
open IN,"<","hzau.txt";
while (<IN>)
{
    chomp;
    if ($_=~/"(http:\/\/[0-9a-z\/\-\.]+[^\/])\/*"/i)
    {
        $hash{$1}=1;
    }
}
foreach $key(keys %hash)
{
    print("$key\n");
}
close IN;
#5
$t++;
print("-"x20,"$t","-"x20,"\n");
@dna=();
while (<FILE1>)
{
    chomp;
    if ($_=~/^\w+/i)
    {
        push(@dna,$_);
    }
}
t5_1(@dna);
t5_2(@dna);
close FILE1;
#6
$t++;
print("-"x20,"$t","-"x20,"\n");
@data=();
$count=0;
while (<FILE2>)
{
    chomp;
    if ($_=~/ATOM\s+\d+\s+CA\s+(\w{3})\s+D\s+(\d{3})\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+).+C/)
    {
        $count++;
        push(@data,[$1,$2,$3,$4,$5]);
    }
}
t6(@data);
#7
$t++;
print("-"x20,"$t","-"x20,"\n");
$a=('a');
$b=('b');
$c=('c');
=pod
void move(int n,char a,char b,char c)
{
    if(n==1)
        printf("\t%c->%c\n",a,c);    //当n只有1个的时候直接从a移动到c
    else
    {
        move(n-1,a,c,b);            //第n-1个要从a通过c移动到b
        printf("\t%c->%c\n",a,c);
        move(n-1,b,a,c);            //n-1个移动过来之后b变开始盘，b通过a移动到c，这边很难理解
    }
}
 
main()
{
    int n;
    printf("请输入要移动的块数：");
    scanf("%d",&n);
    move(n,'a','b','c');
}
=cut
print("请输入要移动的块数:");
chomp($n=<STDIN>);
&t7(($n,$a,$b,$c));

sub t7
{
    if ($_[0]==1)
    {
        print("$_[1]->$_[3]\n");
    }
    else
    {
        t7(($_[0]-1,$_[1],$_[2],$_[3]));
        print("$_[1]->$_[3]\n");
        t7(($_[0]-1,$_[2],$_[1],$_[3]));
    }
}

__FILE1__
>gi|55380579|gb|AE014297.2| Drosophila melanogaster chromosome 3R, complete sequence
GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATGCAATGGTCAGGTAGCGTTGTTCTAAACTCAAGATTT
AAAGGTGAATAGTCCTGTAAGCCCTATAAACATATGTACATAGGTAGGCCAGTACTTAGTACTGGCACAT
GCCGCTGATCTGTTAGTAGATTATCCATTTCCCTTCAGCGCCTACCTGCGTCACCAATGATGAGGTCGAG
ACAGAATCCTACTAGTACCTGCCTCGAGTCGATCGGGCAGAGAGCGAGAAATGGTAAGCAGGTGAGTGAG
CGCAGAGAGCGTCTTTCGACGACTCTTTCGTCGCGAGCAAACAACAAGTAGACGTCGCTCAGACACTGTC
GGCCAGATTCATTTTCCAGAAAGACGTCGTCGCGTTGACAAGCTTAAATTCGTAGCGGGCGCCAGTAGGA
CGACCCAGTGGATATCGTCAGTTGAACCAGGGGAAACGTAGCAGCCCAGTTACATTGCTCGGGAGGGGTA
__FILE2__
ATOM   4542  CA  LYS D 130      27.728  42.781 107.352  1.00 53.17           C  
ATOM   4551  CA  VAL D 131      28.783  40.337 110.073  1.00 46.40           C  
ATOM   4558  CA  LYS D 132      30.084  42.479 112.907  1.00 52.26           C  
ATOM   4567  CA  ALA D 133      31.868  41.767 116.183  1.00 53.70           C  
ATOM   4572  CA  LYS D 134      28.704  41.872 118.303  1.00 71.02           C  
ATOM   4581  CA  SER D 135      26.883  39.124 116.249  1.00 78.80           C  
ATOM   4587  CA  ILE D 136      29.774  36.647 116.820  1.00 93.90           C  
ATOM   4595  CA  VAL D 137      28.482  35.979 120.344  1.00117.21           C  
ATOM   4602  CA  PHE D 138      27.994  33.911 123.464  1.00133.01           C  
ATOM   4613  CA  HIS D 139      26.718  30.262 123.566  1.00139.94           C  
ATOM   4623  CA  ARG D 140      29.353  28.808 121.207  1.00139.51           C  
ATOM   4634  CA  LYS D 141      30.584  25.201 121.443  1.00131.58           C  
ATOM   4642  CA  LYS D 142      33.749  25.778 123.504  1.00113.11           C  
ATOM   4651  CA  ASN D 143      37.231  24.833 122.254  1.00 94.20           C  
ATOM   4659  CA  LEU D 144      36.495  27.447 119.636  1.00 69.72           C  
ATOM   4667  CA  GLN D 145      38.015  30.917 119.225  1.00 40.46           C  
ATOM   4676  CA  TYR D 146      37.259  33.717 116.776  1.00 28.66           C  
ATOM   4688  CA  TYR D 147      39.566  36.141 114.844  1.00 24.27           C  
ATOM   4700  CA  ASP D 148      38.913  38.812 112.220  1.00 32.21           C  
ATOM   4708  CA  ILE D 149      41.467  38.220 109.408  1.00 35.28           C  
ATOM   4716  CA  SER D 150      42.088  39.987 106.101  1.00 32.96           C  
ATOM   4722  CA  ALA D 151      43.967  37.479 103.880  1.00 32.91           C  
ATOM   4727  CA  LYS D 152      44.860  40.086 101.166  1.00 35.10           C  
ATOM   4736  CA  SER D 153      46.596  42.467 103.623  1.00 32.49           C  
ATOM   4742  CA  ASN D 154      47.988  39.768 105.958  1.00 28.75           C  