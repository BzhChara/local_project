use utf8;
use open ":encoding(gbk)",":std";
use List::Util qw/max min sum maxstr minstr shuffle/;

open IN_1,"<","ecoli-soap.contig";
open IN_2,"<","ecoli-soap.scafSeq";

print("contig:\n");
@data=();
while (<IN_1>)
{
    chomp;
    if ($_=~/\>\d+\s+length\s+(\d+)/)
    {
        if ($1>=200)
        {
            push(@data,$1);
        }
    }
}

@DATA=sort{$b<=>$a}(@data);
#print("$DATA[0]\n");
#4660

$sum=0;
foreach $n(0..$#data)
{
    $sum=$sum+$data[$n];
}
#print("$sum\n");
#16349478

$N50=$sum*0.5;
$N90=$sum*0.9;

$count=$n=0;
foreach $n(0..$#DATA)
{
    $count=$count+$DATA[$n];
    if ($count>=$N50)
    {
        printf("N50:$DATA[$n] L50:%d\n",$n+1);
        last;
    }
}

$count=$n=0;
foreach $n(0..$#DATA)
{
    $count=$count+$DATA[$n];
    if ($count>=$N90)
    {
        printf("N90:$DATA[$n] L90:%d\n",$n+1);
        last;
    }
}

$length=@DATA;
print("总碱基数:$sum 总序列数目:$length 最长序列的长度:$DATA[0]\n");

print("scafSeq:\n");
@data=();
$count=0;
while (<IN_2>)
{
    chomp;
    if ($_=~/^\w+/)
    {
        $count=$count+length($_);
    }
    else
    {
        if ($count>=200)
        {
            push(@data,$count);
        }
        $count=0;
    }
}

if ($count>=200)
{
    push(@data,$count);
}
$count=0;

@DATA=sort{$b<=>$a}(@data);

$sum=0;
foreach $n(0..$#data)
{
    $sum=$sum+$data[$n];
}

$N50=$sum*0.5;
$N90=$sum*0.9;

$count=$n=0;
foreach $n(0..$#DATA)
{
    $count=$count+$DATA[$n];
    if ($count>=$N50)
    {
        printf("N50:$DATA[$n] L50:%d\n",$n+1);
        last;
    }
}

$count=$n=0;
foreach $n(0..$#DATA)
{
    $count=$count+$DATA[$n];
    if ($count>=$N90)
    {
        printf("N90:$DATA[$n] L90:%d\n",$n+1);
        last;
    }
}

$length=@DATA;
print("总碱基数:$sum 总序列数目:$length 最长序列的长度:$DATA[0]\n");

close IN_1;
close IN_2;
