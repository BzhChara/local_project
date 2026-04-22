use utf8;
use open ":encoding(gbk)",":std";

#查看数组内容
sub test
{
    foreach $i(0..$#_)
    {
        print("@{$_[$i]}\n");
    }
}

open IN,"<","2a07.pdb";
open OUT_F,">","2A07_F.pdb";
open OUT_F_DNA,">","2A07_F_DNA.pdb";
while (<IN>)
{
    chomp;
    if ($_=~/(ATOM\s+\d+\s+[\w\d]{1,3}\s+\w{3}\s+F.+\w)/)
    {
        print{OUT_F}("$1\n");
    }
    if ($_=~/(ATOM\s+\d+\s+[\w\d\']{1,3}\s+D[A|G|C|T].+\w)/)
    {
        print{OUT_F_DNA}("$1\n");
    }
}
close IN;
close OUT_F;
close OUT_F_DNA;

$tmp=();
%hash=('A'=>106,'C'=>135,'D'=>163,'E'=>194,'F'=>197,'G'=>84,'H'=>184,'I'=>169,'K'=>205,'L'=>164,'M'=>188,'N'=>157,'P'=>136,'Q'=>198,'R'=>248,'S'=>130,'T'=>142,'V'=>142,'W'=>227,'Y'=>222);
open IN,"<","2A07_F.dssp";
open OUT,">","2A07_F_rasa.txt";
while (<IN>)
{
    chomp;
    if ($_=~/\s+\d+\s+(\d{3})\s+F\s+(\w).{20}\s+(\d+)/)
    {
        $tmp=$3/$hash{$2};
        print{OUT}("$1  $2  $tmp\n");
    }
}
close IN;
close OUT;

open IN_1,"<","2A07_F_rasa.txt";
open IN_2,"<","2A07_F.pdb";
open IN_3,"<","2A07_F_DNA.pdb";
open OUT,">","2A07_F_binding.txt";
@data_1=();
@data_2=();
@data_3=();
while (<IN_1>)
{
    chomp;
    push(@data_1,[split('\s+',$_)]);
}
#test(@data_1);

while (<IN_2>)
{
    chomp;
    push(@data_2,[split('\s+',$_)]);
}
#test(@data_2);

while (<IN_3>)
{
    chomp;
    push(@data_3,[split('\s+',$_)]);
}
#test(@data_3);

while ($x<=$#data_1)
{
    if ($data_1[$x][2]>0.1)
    {
        foreach $y(0..$#data_2)
        {
            if ($data_1[$x][0]==$data_2[$y][5])
            {
                foreach $z(0..$#data_3)
                {
                    $length=sqrt((($data_2[$y][6]-$data_3[$z][6])**2)+(($data_2[$y][7]-$data_3[$z][7])**2)+(($data_2[$y][8]-$data_3[$z][8])**2));
                    if ($length<4.5)
                    {
                        print{OUT}("1\t$data_1[$x][0]\t$data_1[$x][1]\n");
                        goto part1;
                    }
                }
            }
        }
        print{OUT}("0\t$data_1[$x][0]\t$data_1[$x][1]\n");
    }
    else
    {
        print{OUT}("0\t$data_1[$x][0]\t$data_1[$x][1]\n");
    }
    part1:$x++;
}
close IN_1;
close IN_2;
close IN_3;
close OUT;