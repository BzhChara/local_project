#!/usr/bin/perl
use Inline::Files;
use utf8;
use open ":encoding(gbk)",":std";

#不能读入FILE请安装cpanm install Inline::Files
#如果不能读入数据，请将文件改为不带中文的名称

$t=0;#题号

#1
%hash=();$t++;
print("*"x20,"$t","*"x20,"\n");
while (<FILE1>)
{
    chomp;
    $_=~/(\w*)\s*=>\s*(\w*)/;
    $hash{$1}=$2;
}
print("-"x5,"(1)","-"x5,"\n");
foreach $key(sort{$b cmp $a} keys %hash)
{
    print("$key=>$hash{$key}\n");
}
print("-"x5,"(2)","-"x5,"\n");
foreach $key(sort{$hash{$b} cmp $hash{$a}} keys %hash)
{
    print("$key=>$hash{$key}\n");
}
close FILE1;
#2
%hash=();$t++;
print("*"x20,"$t","*"x20,"\n");
while (<FILE2>)
{
    chomp;
    $_=~/(\d*)\s*=>\s*(\d*)/;
    $hash{$1}=$2;
}
foreach $key(sort{($hash{$a} <=> $hash{$b}) || ($b <=> $a)} keys %hash)
{
    print("$key=>$hash{$key}\n");
}
close FILE2;
#3
%hash=();$t++;
print("*"x20,"$t","*"x20,"\n");
while (<FILE3>)
{
    chomp;
    $_=~/(\d*)\s*(\d*)\s*(\d*)/;
    #print("$1-$2-$3\n");
    $hash{$1}{$2}=$3;
}
@p=(),%q=();
foreach $k1(keys %hash)
{
    push(@p,$k1);
	foreach $k2(keys %{$hash{$k1}}) 
    {
		push(@p,$k2);
	}
}
foreach (@p)
{
    $q{$_}=1;
}
@portein=sort{$a<=>$b}keys %q;
$length=@portein;
@matrix=();
push(@matrix,["000",@portein]);
foreach (0..$length-1)
{
    push(@matrix,[$portein[$_],"000","000","000","000","000","000","000","000","000","000","000","000","000","000"]);
}
foreach $k1(keys %hash)
{
	foreach $k2(keys %{$hash{$k1}}) 
    {
        foreach $i(1..$length)
        {
            foreach $j(1..$length)
            {
                if ($matrix[$i][0]==$k1 && $matrix[0][$j]==$k2)
                {
                    $matrix[$i][$j]=$hash{$k1}{$k2};
                }
                if ($matrix[$i][0]==$k2 && $matrix[0][$j]==$k1)
                {
                    $matrix[$i][$j]=$hash{$k1}{$k2};
                }
            }
        }
	}
}
foreach (0..$#matrix)
{
    print("@{$matrix[$_]}\n");
}
close FILE3;
#4
%hash=();$t++;
print("*"x20,"$t","*"x20,"\n");
print("-"x5,"(a)","-"x5,"\n");
@rna=();
@total=();
$num=0;
while (<FILE4>)
{
    chomp;
    $_=~/>\d*_*\d*_*([0-9]*)/;
    push(@total,$1);
    $_=~/(\w*)/;
    push(@rna,split('',$1));
    if ($_=~/>/)
    {
        $num++;
    }
}
$n=0;
foreach (0..$#rna)
{
    print("$rna[$_]");
    $n++;
    if ($n==80)
    {
        print("\n");
        $n=0;
    }
}
print("\n");
print("-"x5,"(b)","-"x5,"\n");
$count=0;
foreach (0..$#rna-1)
{
    if ($rna[$_] eq 'G' || $rna[$_] eq 'C')
    {
        $count++;
    }
}
printf("%f%\n",$count/(scalar(@rna)-1)*100);
print("-"x5,"(c)","-"x5,"\n");
$count=0;
print("$num\n");
print("-"x5,"(d)","-"x5,"\n");
$n=0;
foreach (0..$#total)
{
    $n=$n+$total[$_];
}
print("$n\n");
close FILE4;
#5
%hash=();$t++;
print("*"x20,"$t","*"x20,"\n");
print("-"x5,"(a)","-"x5,"\n");
@Query=();
@Subject=();
while (<FILE5>)
{
    chomp;
    if ($_=~/Query\s*\d*\s*([A-Z,-]*)\s*\d*.*/)
    {
        push(@Query,split('',$1));
    }
    if ($_=~/Sbjct\s*\d*\s*([*,A-Z,-]*)\s*\d*/)
    {
        push(@Subject,split('',$1));
    }
}
print(@Query,"\n");
print(@Subject,"\n");
print("-"x5,"(b)","-"x5,"\n");
$query=join('',@Query);
$subject=join('',@Subject);
$query=~tr/-//d;
$subject=~tr/-//d;
print("$query\n");
print("$subject\n");
print("-"x5,"(c)","-"x5,"\n");
$length_1=length($query);
print("Query:$length_1\n");
$length_2=length($subject);
print("Subject:$length_2\n");
$count_1=$count_2=0;
@type=qw/AA AG AC AT GA GG GC GT CA CG CC CT TA TG TC TT/;
foreach $i(0..$length_1-1)
{
    foreach $j(0..$#type-1)
    {
        if (substr($query,$i,2) eq $type[$j])
        {
            $count_1++;
        }
    }
}
foreach $i(0..$length_2-1)
{
    foreach $j(0..$#type-1)
    {
        if (substr($subject,$i,2) eq $type[$j])
        {
            $count_2++;
        }
    }
}
printf("Query:%f%\n",$count_1/$length_1*100);
printf("Subject:%f%\n",$count_2/$length_2*100);
close FILE5;
#6
%hash=();$t++;
print("*"x20,"$t","*"x20,"\n");
print("类别 种类\n");
while (<FILE6>)
{
    chomp;
    ($who,$rest)=split('/\s*/',$_,2);
    @fileds=split('/\s*/',$rest);
    $HoA{$who}=[@fileds];
}
foreach $family(keys %HoA)
{
    print("$family @{$HoA{$family}}\n");
}
close FILE6;
#7
%hash=();$t++;
print("*"x20,"$t","*"x20,"\n");
open IN,"<","bicluster.txt";
open OUT,">","bicluster_output.txt";
@bicluster=();
while (<IN>)
{
    chomp;
    push(@bicluster,split(' ',$_));
}
foreach (0..$#bicluster)
{
    if ($bicluster[$_]=~/(bicluster\d*)/)
    {
        if ($_>0)
        {
            print{OUT}("\n\n");
        }
        print{OUT}("$1 ");
    }
    if ($bicluster[$_]=~/(.{6})[(]\d*[,]\d*[,]\d*[,]\d*[,](.{1,10})[)]/)
    {
        print{OUT}("$1:$2 ");
    }
}
print("详细请见文件\n");
close IN;
close OUT;
#8
%hash=();$t++;
print("*"x20,"$t","*"x20,"\n");
open IN,"<","anjisuan.txt";
open OUT,">","nonredundant.txt";
@portein=();
while (<IN>)
{
    chomp;
    if ($_=~/(\w+)/)
    {
        push(@portein,$1);
    }
}
for ($i=0;$i<$#portein-1;$i=$i+2)
{
    #print("$portein[$i].$portein[$i+1]\n");
    $hash{join('',($portein[$i],$portein[$i+1]))}+=1;
}
foreach $key(sort{$hash{$a}<=>$hash{$b}} keys %hash)
{
    print{OUT}("$key $hash{$key}\n");
}
print("详细请见文件\n");
close IN;
close OUT;
open IN,"<","nonredundant.txt";
%hash=%num=();
$min=$max=$count_1=$count_2=0;
print("请输入上限:\n");
chomp($max=<STDIN>);
print("请输入下限\n");
chomp($min=<STDIN>);
while (<IN>)
{
    chomp;
    if ($_=~/(\w+)\s(\d+)/)
    {
        $hash{$1}=$2;
    }
}
foreach $key(sort{$hash{$a}<=>$hash{$b}} keys %hash)
{
    if ($hash{$key}>=$min && $hash{$key}<=$max)
    {
        $count_1=$count_1+length($key)*$hash{$key};
        @data=split('',$key);
        foreach (0..$#data)
        {
            $num{$data[$_]}+=1*$hash{$key};
        }
    }
}
foreach $key(keys %num)
{
    print("$key $num{$key} ");
    printf("%f%\n",$num{$key}/$count_1*100);
}
close IN;

__FILE1__
Donald => Knuth
Alan   => Turing
John   => Neumann
Sum    => Pearsom
Soddy  => Toleis
__FILE2__
12  =>  50;
90  =>  50;
60  =>  7;
49  =>  100;
8   =>  100;
__FILE3__
123 345 3
233 333 23
455	344 45
344 444 65
444 449 91
421 457 34
221 409 45
467 222 40
__FILE4__
>2_24_3
AAAAAACAACCTCTCTACCTGTTC
>7_17_6
AAAAAACAAGTAGATCA
>8_24_1
AAAAAACAATTAACTGTGGACGGA
>9_18_2
AAAAAACACAATCAAATA
>11_24_2
AAAAAACAGACTGCAGTTGACGAT
>12_24_62
AAAAAACAGGCTGAGACGACGGAA
>13_25_6
AAAAAACAGGCTGAGACGACGGAAC
>14_24_7
AAAAAACAGGCTGAGACGATGGAA
>15_24_11
AAAAAACAGGCTGAGATGACGGAA
>16_24_31
AAAAAACAGGGCTGAGACGACGGA
>18_15_5
AAAAAACATTTTCCT
>19_24_6
AAAAAACCACGGACCAAGAGGAGC
>20_24_2
AAAAAACCACGGACTAAGAGAAGC
>21_24_6
AAAAAACCACGGACTAAGAGGAGC
>22_24_4
AAAAAACCACGGGCCAAGAGAAGC
>23_24_18
AAAAAACCACGGGCCAAGAGGAGC
>24_24_2
AAAAAACCATAACCTACCTCCACC
>25_22_1
AAAAAACCATCTTTCAATTTCT
>27_17_2
AAAAAACCCATTTTAAA
>30_24_3
AAAAAACCGATAGTTGCGAATTCA
>33_23_5
AAAAAACGGAAAGAGGTAGTCAA
>34_24_5
AAAAAACGGACTTCAAAGTAGATT
>35_24_7
AAAAAACGGGCTGAGACGACGGAA
>36_25_4
AAAAAACGGGCTGAGACGACGGAAC
>39_24_14
AAAAAACGTTGGGCACAGAAGATA
>40_24_5
AAAAAACTAAATGCAAATGACGGC
>41_15_4
AAAAAACTAAGAAAA
>44_23_1
AAAAAACTACAGTAAATCACCAT
>45_15_19
AAAAAACTACTGAAT
>46_24_11
AAAAAACTACTGGCCAAGAGGAGC
>48_24_7
AAAAAACTATCGGATCTCCTCCAT
>49_24_7
AAAAAACTATGAACAACTGAACGC
>51_15_1
AAAAAACTCCTATTT
>52_23_8
AAAAAACTCGAGACCTCTGATTA
>54_23_3
AAAAAACTCTACAACAGATGATG
>55_22_9
AAAAAACTCTGATCAAGCATGA
>57_25_2
AAAAAACTCTGATCAAGCATGACGA
>58_24_22
AAAAAACTCTGTTCACGGTATATA
>60_24_23
AAAAAACTGAACAAAGCGGAAGAC
>61_24_2
AAAAAACTGAACCGACTGATAAGC
>63_23_6
AAAAAACTGAAGGAGTCTGTCAA
>65_24_4
AAAAAACTGGCTGAGACGACGGAA
>66_25_5
AAAAAACTGGCTGAGACGACGGAAC
>67_24_2
AAAAAACTGGCTGAGACGACGGAT
>68_23_7
AAAAAACTGTAGCAAGCGAGAAG
>69_24_2
AAAAAACTGTATAATTCGCTATAC
>72_23_1
AAAAAACTTCCTTATCTTAACAA
>74_24_5
AAAAAACTTGGACAAGGCGGCTCA
>75_24_11
AAAAAACTTGTAGAAAGCAGAAGA
>78_24_6
AAAAAACTTTTAACAGATGGCATA
__FILE5__
Query  333     WISLGAVKRLQEQNAIDMEIIVNPKT------------------------------LDGG  362
               W++L A+K+L E +A+ MEII NPK                               +DG 
Sbjct  751020  WVNLKAIKKLVEADALKMEIIPNPKVPAWLMHSFLFLQFFHSLLSNQRKTANSIQEVDG-  750844

Query  363     LNVIQLETAVGAAIKS--------------------------------FENSLGINVPRS  390
               + V+QLETA GAAI+                                 F+N++G+NVPRS
Sbjct  750843  VKVLQLETAAGAAIRVYMPVLVARFEI*NANSNLFLTPFSCDSVSKQFFDNAIGVNVPRS  750664

Query  391     RFLPVKTTSDLLLV-------------------------------MSNLYSLNAGSLTMS  419
               RFLPVK +SDLLLV                                S+LY+L  G +T +
Sbjct  750663  RFLPVKASSDLLLVQVNSIHHV*QRNRVTNIF*NI*LMI*VFFLPQSDLYTLVDGFVTRN  750484

Query  420     EKREFPTVPLVKLGSSFTKVQD-----------------------------YLRRFESIP  450
               + R  P+ P ++LG  F KV +                             +L RF+SIP
Sbjct  750483  KARTNPSNPSIELGPEFKKVNEENDFLTSPICHFYLDS*ASLSVSV*QVATFLSRFKSIP  750304

Query  451     DMLELDHLTVSGDVTFGKNVSLK  473
                ++ELD L VSGDV FG ++ LK
Sbjct  750303  SIVELDSLKVSGDVWFGSSIVLK  750235
__FILE6__
plant	arabidopsis oryza rape maize
animal	swine cattle
human	you me her