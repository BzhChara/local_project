use utf8;
use open ":encoding(gbk)",":std";
use List::Util qw/max min sum maxstr minstr shuffle/;

#序列（长度为15）
@DNA_1=qw/A A C G T A C T C A A G T C T/;
@DNA_2=qw/T C G T A C T C T A A C G A T/;

print("-"x20,"Needleman-Wunsch算法:全局比对","-"x20,"\n");

#评分规则1
$match=9;
$mismatch=-6;
$insertion=$deletion=-2;

#构建打分矩阵
@matrix=();
push(@matrix,[' ',' ',@DNA_1]);
push(@matrix,[' ',0,-2,-4,-6,-8,-10,-12,-14,-16,-18,-20,-22,-24,-26,-28,-30]);
$n=0;
foreach (0..$#DNA_2)
{
    $n=$n-2;
    push(@matrix,[$DNA_2[$_],$n]);
}

#评分
@tmp=(0,0,0);
for ($n=2;$n<17;$n++)
{
    foreach $m(2..16)
    {
        $tmp[0]=$matrix[$m-1][$n]-2;
        $tmp[1]=$matrix[$m][$n-1]-2;
        if ($matrix[$m][0] eq $matrix[0][$n])
        {
            $tmp[2]=$matrix[$m-1][$n-1]+9;
        }
        else
        {
            $tmp[2]=$matrix[$m-1][$n-1]-6;
        }
        $matrix[$m][$n]=max @tmp;
        $tmp[0]=$matrix[$n-1][$m]-2;
        $tmp[1]=$matrix[$n][$m-1]-2;
        if ($matrix[0][$m] eq $matrix[$n][0])
        {
            $tmp[2]=$matrix[$n-1][$m-1]+9;
        }
        else
        {
            $tmp[2]=$matrix[$n-1][$m-1]-6;
        }
        $matrix[$n][$m]=max @tmp;
    }
}

#打印
foreach $i(0..16)
{
    foreach $j(0..16)
    {
        print("$matrix[$i][$j]");
        if ($j==16)
        {
            print("\n");
        }
        if ($j!=16)
        {
            print("\t");
        }
    }
}

#构建回溯函数
sub BZH
{
    open OUT,">","Needleman-Wunsch.txt";
    my ($i,$j,@data)=@_;
    my @tmp=();
    if (($data[$i-1][$j]-2) == $data[$i][$j])
    {
        push(@tmp,[($i-1),$j]);
        #printf("%d:%d",$i-1,$j);
    }
    if (($data[$i][$j-1]-2) == $data[$i][$j])
    {
        push(@tmp,[$i,($j-1)]);
        #printf("%d:%d",$i,$j-1);
    }
    if (($data[$i-1][$j-1]-6)==$data[$i][$j])
    {
        if ($data[$i][0] ne $data[0][$j])
        {
            push(@tmp,[($i-1),($j-1)]);
            #printf("%d:%d",$i-1,$j-1);
        }
    }
    if (($data[$i-1][$j-1]+9)==$data[$i][$j])
    {
        if ($data[$i][0] eq $data[0][$j])
        {
            push(@tmp,[($i-1),($j-1)]);
            #printf("%d:%d",$i-1,$j-1);
        }
    }
    foreach my $n(0..$#tmp)
    {
        print{OUT}("$tmp[$n][0]:$tmp[$n][1]");
        print{OUT}(" "x10);
    }
    print{OUT}("\n");
    close OUT;

part1:open IN,"<","Needleman-Wunsch.txt";
    my @path=();
    my @array=<IN>;
    my $tail=$array[$#array];
    while ($tail=~/(\d+)\:(\d+)/g)
    {
        chomp;
        push(@path,[$1,$2]);
    }
    close IN;
    foreach $m(0..$#path)
    {
        open OUT,">>","Needleman-Wunsch.txt";
        @tmp=();
        $i=$path[$m][0];
        $j=$path[$m][1];
        if (($data[$i-1][$j]-2) == $data[$i][$j])
        {
            push(@tmp,[($i-1),$j]);
            #printf("%d:%d",$i-1,$j);
        }
        if (($data[$i][$j-1]-2) == $data[$i][$j])
        {
            push(@tmp,[$i,($j-1)]);
            #printf("%d:%d",$i,$j-1);
        }
        if (($data[$i-1][$j-1]-6)==$data[$i][$j])
        {
            if ($data[$i][0] ne $data[0][$j])
            {
                push(@tmp,[($i-1),($j-1)]);
                #printf("%d:%d",$i-1,$j-1);
            }
        }
        if (($data[$i-1][$j-1]+9)==$data[$i][$j])
        {
            if ($data[$i][0] eq $data[0][$j])
            {
                push(@tmp,[($i-1),($j-1)]);
                #printf("%d:%d",$i-1,$j-1);
            }
        }
        foreach $n(0..$#tmp)
        {
            print{OUT}("$tmp[$n][0]:$tmp[$n][1]");
            print{OUT}(" "x10);
        }
    }
    print{OUT}("\n");
    close OUT;
    if ($i>=2 && $j>=2)
    {
        goto part1;
    }
    else
    {
        return;
    }
}

BZH(16,16,@matrix);

#去重
open IN,"<","Needleman-Wunsch.txt";
open OUT,">","NewNeedleman-Wunsch.txt";
%hash=();
print{OUT}("16:16\n");
while (<IN>)
{
    chomp;
    while ($_=~/(\d+\:\d+)/g)
    {
        $hash{$1}=1;
    }
    foreach $key(keys %hash)
    {
        print{OUT}("$key\t");
    }
    print{OUT}("\n");
    %hash=();
}
close IN;
close OUT;

=pod
open IN,"<","New.txt";
while (<IN>)
{
    chomp;
    while ($_=~/(\d+)\:(\d+)/g)
    {
        if ($1!=0 && $2!=0)
        {
            $matrix[$1][$2]='XXX';
        }
    }
}
close IN;
=cut

#储存匹配结果
open IN,"<","NewNeedleman-Wunsch.txt";
@array=<IN>;
@{$S[0][0]}=();
for ($n=$#array;$n>0;$n--)
{
    chomp;
    @array_1=();
    @array_2=();
    while ($array[$n]=~/(\d+)\:(\d+)/g)
    {
        if ($1!=0 && $2!=0)
        {
            push(@array_1,[$1,$2]);
        }
    }
    while ($array[$n-1]=~/(\d+)\:(\d+)/g)
    {
        push(@array_2,[$1,$2]);
    }
    foreach $i(0..$#array_1)
    {
        foreach $j(0..$#array_2)
        {
            if (($array_2[$j][0]-$array_1[$i][0])==1 && ($array_2[$j][1]-$array_1[$i][1])==1 && ($matrix[$array_1[$i][0]][$array_1[$i][1]]+9)==($matrix[$array_2[$j][0]][$array_2[$j][1]]) && $matrix[$array_2[$j][0]][0] eq $matrix[0][$array_2[$j][1]])
            {
                @{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}=$DNA_2[$array_1[$i][0]-1].$DNA_1[$array_1[$i][1]-1];
                #print("@{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}\n");
            }
            if (($array_2[$j][0]-$array_1[$i][0])==1 && ($array_2[$j][1]-$array_1[$i][1])==1 && ($matrix[$array_1[$i][0]][$array_1[$i][1]]-6)==($matrix[$array_2[$j][0]][$array_2[$j][1]]) && $matrix[$array_2[$j][0]][0] ne $matrix[0][$array_2[$j][1]])
            {
                @{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}=$DNA_2[$array_1[$i][0]-1].$DNA_1[$array_1[$i][1]-1];
                #print("@{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}\n");
            }
            if (($array_2[$j][0]-$array_1[$i][0])==1 && ($array_2[$j][1]-$array_1[$i][1])==0 && ($matrix[$array_1[$i][0]][$array_1[$i][1]]-2)==($matrix[$array_2[$j][0]][$array_2[$j][1]]))
            {
                @{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}=$DNA_2[$array_1[$i][0]-1].'-';
                #print("@{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}\n");
            }
            if (($array_2[$j][0]-$array_1[$i][0])==0 && ($array_2[$j][1]-$array_1[$i][1])==1 && ($matrix[$array_1[$i][0]][$array_1[$i][1]]-2)==($matrix[$array_2[$j][0]][$array_2[$j][1]]))
            {
                @{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}='-'.$DNA_1[$array_1[$i][1]-1];
                #print("@{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}\n");
            }
        }
    }
}
close IN;

#结果输出(还没有想到好的输出方法,所以方法很蠢)
open IN,"<","NewNeedleman-Wunsch.txt";
@array=<IN>;
@array_1=();
@array_2=();
@array_4=();
@array_5=();
@array_6=();
@array_7=();
@array_8=();
@array_9=();
@array_10=();
@array_11=();
@array_12=();
@array_13=();
@array_14=();
@array_15=();
@array_16=();
@array_17=();
@array_18=();
@array_19=();
@array_20=();
$count=0;
while ($array[$#array]=~/(\d+)\:(\d+)/g)
{
    if ($1!=0 && $2!=0)
    {
        push(@array_1,[$1,$2]);
    }
}
while ($array[$#array-1]=~/(\d+)\:(\d+)/g)
{
    push(@array_2,[$1,$2]);
}
while ($array[$#array-2]=~/(\d+)\:(\d+)/g)
{
    push(@array_3,[$1,$2]);
}
while ($array[$#array-3]=~/(\d+)\:(\d+)/g)
{
    push(@array_4,[$1,$2]);
}
while ($array[$#array-4]=~/(\d+)\:(\d+)/g)
{
    push(@array_5,[$1,$2]);
}
while ($array[$#array-5]=~/(\d+)\:(\d+)/g)
{
    push(@array_6,[$1,$2]);
}
while ($array[$#array-6]=~/(\d+)\:(\d+)/g)
{
    push(@array_7,[$1,$2]);
}
while ($array[$#array-7]=~/(\d+)\:(\d+)/g)
{
    push(@array_8,[$1,$2]);
}
while ($array[$#array-8]=~/(\d+)\:(\d+)/g)
{
    push(@array_9,[$1,$2]);
}
while ($array[$#array-9]=~/(\d+)\:(\d+)/g)
{
    push(@array_10,[$1,$2]);
}
while ($array[$#array-10]=~/(\d+)\:(\d+)/g)
{
    push(@array_11,[$1,$2]);
}
while ($array[$#array-11]=~/(\d+)\:(\d+)/g)
{
    push(@array_12,[$1,$2]);
}
while ($array[$#array-12]=~/(\d+)\:(\d+)/g)
{
    push(@array_13,[$1,$2]);
}
while ($array[$#array-13]=~/(\d+)\:(\d+)/g)
{
    push(@array_14,[$1,$2]);
}
while ($array[$#array-14]=~/(\d+)\:(\d+)/g)
{
    push(@array_15,[$1,$2]);
}
while ($array[$#array-15]=~/(\d+)\:(\d+)/g)
{
    push(@array_16,[$1,$2]);
}
while ($array[$#array-16]=~/(\d+)\:(\d+)/g)
{
    push(@array_17,[$1,$2]);
}
while ($array[$#array-17]=~/(\d+)\:(\d+)/g)
{
    push(@array_18,[$1,$2]);
}
while ($array[$#array-18]=~/(\d+)\:(\d+)/g)
{
    push(@array_19,[$1,$2]);
}
while ($array[$#array-19]=~/(\d+)\:(\d+)/g)
{
    push(@array_20,[$1,$2]);
}
foreach $i(0..$#array_1)
{
    foreach $j(0..$#array_2)
    {
        foreach $m(0..$#array_3)
        {
            foreach $m1(0..$#array_4)
            {
                foreach $m2(0..$#array_5)
                {
                    foreach $m3(0..$#array_6)
                    {
                        foreach $m4(0..$#array_7)
                        {
                            foreach $m5(0..$#array_8)
                            {
                                foreach $m6(0..$#array_9)
                                {
                                    foreach $m7(0..$#array_10)
                                    {
                                        foreach $m8(0..$#array_11)
                                        {
                                            foreach $m9(0..$#array_12)
                                            {
                                                foreach $m10(0..$#array_13)
                                                {
                                                    foreach $m11(0..$#array_14)
                                                    {
                                                        foreach $m12(0..$#array_15)
                                                        {
                                                            foreach $m13(0..$#array_16)
                                                            {
                                                                foreach $m14(0..$#array_17)
                                                                {
                                                                    foreach $m15(0..$#array_18)
                                                                    {
                                                                        foreach $m16(0..$#array_19)
                                                                        {
                                                                            foreach $m17(0..$#array_20)
                                                                            {
                                                                                if (@{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]})
                                                                                {       
                                                                                    if (@{$S[$array_2[$j][0].$array_2[$j][1]][$array_3[$m][0].$array_3[$m][1]]})
                                                                                    {
                                                                                        if (@{$S[$array_3[$m][0].$array_3[$m][1]][$array_4[$m1][0].$array_4[$m1][1]]})
                                                                                        {
                                                                                            if (@{$S[$array_4[$m1][0].$array_4[$m1][1]][$array_5[$m2][0].$array_5[$m2][1]]})
                                                                                            {
                                                                                                if (@{$S[$array_5[$m2][0].$array_5[$m2][1]][$array_6[$m3][0].$array_6[$m3][1]]})
                                                                                                {
                                                                                                    if (@{$S[$array_6[$m3][0].$array_6[$m3][1]][$array_7[$m4][0].$array_7[$m4][1]]})
                                                                                                    {
                                                                                                        if (@{$S[$array_7[$m4][0].$array_7[$m4][1]][$array_8[$m5][0].$array_8[$m5][1]]})
                                                                                                        {
                                                                                                            if (@{$S[$array_8[$m5][0].$array_8[$m5][1]][$array_9[$m6][0].$array_9[$m6][1]]})
                                                                                                            {
                                                                                                                if (@{$S[$array_9[$m6][0].$array_9[$m6][1]][$array_10[$m7][0].$array_10[$m7][1]]})
                                                                                                                {
                                                                                                                    if (@{$S[$array_10[$m7][0].$array_10[$m7][1]][$array_11[$m8][0].$array_11[$m8][1]]})
                                                                                                                    {
                                                                                                                        if (@{$S[$array_11[$m8][0].$array_11[$m8][1]][$array_12[$m9][0].$array_12[$m9][1]]})
                                                                                                                        {
                                                                                                                            if (@{$S[$array_12[$m9][0].$array_12[$m9][1]][$array_13[$m10][0].$array_13[$m10][1]]})
                                                                                                                            {
                                                                                                                                if (@{$S[$array_13[$m10][0].$array_13[$m10][1]][$array_14[$m11][0].$array_14[$m11][1]]})
                                                                                                                                {
                                                                                                                                    if (@{$S[$array_14[$m11][0].$array_14[$m11][1]][$array_15[$m12][0].$array_15[$m12][1]]})
                                                                                                                                    {
                                                                                                                                        if (@{$S[$array_15[$m12][0].$array_15[$m12][1]][$array_16[$m13][0].$array_16[$m13][1]]})
                                                                                                                                        {
                                                                                                                                            if (@{$S[$array_16[$m13][0].$array_16[$m13][1]][$array_17[$m14][0].$array_17[$m14][1]]})
                                                                                                                                            {
                                                                                                                                                if (@{$S[$array_17[$m14][0].$array_17[$m14][1]][$array_18[$m15][0].$array_18[$m15][1]]})
                                                                                                                                                {
                                                                                                                                                    if (@{$S[$array_18[$m15][0].$array_18[$m15][1]][$array_19[$m16][0].$array_19[$m16][1]]})
                                                                                                                                                    {
                                                                                                                                                        if (@{$S[$array_19[$m16][0].$array_19[$m16][1]][$array_20[$m17][0].$array_20[$m17][1]]})
                                                                                                                                                        {
                                                                                                                                                            $count++;
                                                                                                                                                            print("$count:\n");
                                                                                                                                                            print("@{$S[$array_1[$i][0].$array_1[$i][1]][$array_2[$j][0].$array_2[$j][1]]}  ");
                                                                                                                                                            print("@{$S[$array_2[$j][0].$array_2[$j][1]][$array_3[$m][0].$array_3[$m][1]]}  ");
                                                                                                                                                            print("@{$S[$array_3[$m][0].$array_3[$m][1]][$array_4[$m1][0].$array_4[$m1][1]]}  ");
                                                                                                                                                            print("@{$S[$array_4[$m1][0].$array_4[$m1][1]][$array_5[$m2][0].$array_5[$m2][1]]}  ");
                                                                                                                                                            print("@{$S[$array_5[$m2][0].$array_5[$m2][1]][$array_6[$m3][0].$array_6[$m3][1]]}  ");
                                                                                                                                                            print("@{$S[$array_6[$m3][0].$array_6[$m3][1]][$array_7[$m4][0].$array_7[$m4][1]]}  ");
                                                                                                                                                            print("@{$S[$array_7[$m4][0].$array_7[$m4][1]][$array_8[$m5][0].$array_8[$m5][1]]}  ");
                                                                                                                                                            print("@{$S[$array_8[$m5][0].$array_8[$m5][1]][$array_9[$m6][0].$array_9[$m6][1]]}  ");
                                                                                                                                                            print("@{$S[$array_9[$m6][0].$array_9[$m6][1]][$array_10[$m7][0].$array_10[$m7][1]]}  ");
                                                                                                                                                            print("@{$S[$array_10[$m7][0].$array_10[$m7][1]][$array_11[$m8][0].$array_11[$m8][1]]}  ");
                                                                                                                                                            print("@{$S[$array_11[$m8][0].$array_11[$m8][1]][$array_12[$m9][0].$array_12[$m9][1]]}  ");
                                                                                                                                                            print("@{$S[$array_12[$m9][0].$array_12[$m9][1]][$array_13[$m10][0].$array_13[$m10][1]]}  ");
                                                                                                                                                            print("@{$S[$array_13[$m10][0].$array_13[$m10][1]][$array_14[$m11][0].$array_14[$m11][1]]}  ");
                                                                                                                                                            print("@{$S[$array_14[$m11][0].$array_14[$m11][1]][$array_15[$m12][0].$array_15[$m12][1]]}  ");
                                                                                                                                                            print("@{$S[$array_15[$m12][0].$array_15[$m12][1]][$array_16[$m13][0].$array_16[$m13][1]]}  ");
                                                                                                                                                            print("@{$S[$array_16[$m13][0].$array_16[$m13][1]][$array_17[$m14][0].$array_17[$m14][1]]}  ");
                                                                                                                                                            print("@{$S[$array_17[$m14][0].$array_17[$m14][1]][$array_18[$m15][0].$array_18[$m15][1]]}  ");
                                                                                                                                                            print("@{$S[$array_18[$m15][0].$array_18[$m15][1]][$array_19[$m16][0].$array_19[$m16][1]]}  ");
                                                                                                                                                            print("@{$S[$array_19[$m16][0].$array_19[$m16][1]][$array_20[$m17][0].$array_20[$m17][1]]}  ");
                                                                                                                                                            print("\n");
                                                                                                                                                        }
                                                                                                                                                    }
                                                                                                                                                }
                                                                                                                                            }
                                                                                                                                        }
                                                                                                                                    }
                                                                                                                                }
                                                                                                                            }
                                                                                                                        }
                                                                                                                    }
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
        }
    }
}                                                                             
close IN;

print("-"x20,"Smith-Waterman算法:局部比对","-"x20,"\n");

#评分规则2
$match=9;
$mismatch=-3;
$insertion=$deletion=-2;

#构建打分矩阵
@matrix=();
push(@matrix,[' ',' ',@DNA_1]);
push(@matrix,[' ',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
foreach (0..$#DNA_2)
{
    push(@matrix,[$DNA_2[$_],0]);
}

#评分
@tmp=(0,0,0,0);#多加一个0用于评判正负
for ($n=2;$n<17;$n++)
{
    foreach $m(2..16)
    {
        $tmp[0]=$matrix[$m-1][$n]-2;
        $tmp[1]=$matrix[$m][$n-1]-2;
        if ($matrix[$m][0] eq $matrix[0][$n])
        {
            $tmp[2]=$matrix[$m-1][$n-1]+9;
        }
        else
        {
            $tmp[2]=$matrix[$m-1][$n-1]-6;
        }
        $matrix[$m][$n]=max @tmp;
        $tmp[0]=$matrix[$n-1][$m]-2;
        $tmp[1]=$matrix[$n][$m-1]-2;
        if ($matrix[0][$m] eq $matrix[$n][0])
        {
            $tmp[2]=$matrix[$n-1][$m-1]+9;
        }
        else
        {
            $tmp[2]=$matrix[$n-1][$m-1]-6;
        }
        $matrix[$n][$m]=max @tmp;
    }
}

#打印
foreach $i(0..16)
{
    foreach $j(0..16)
    {
        print("$matrix[$i][$j]");
        if ($j==16)
        {
            print("\n");
        }
        if ($j!=16)
        {
            print("\t");
        }
    }
}

#找到最大值的坐标
$Max=$matrix[1][1];
$row=0;
$col=0;
foreach $i(1..$#matrix)
{
    foreach $j(1..$#{$matrix[$i]})
    {
        if ($matrix[$i][$j]>$Max)
        {
            $Max=$matrix[$i][$j];
            $row=$i;
            $col=$j;
        }
    }
}
#93 16 14
#print("$Max $row $col\n");

#由于局部比对不需要找到全部路径，找到任意一条即可,故回溯函数的构建不需要太严谨，能得出一条即可


open OUT,">","Smith-Waterman.txt";
print{OUT}("16:14\n");
close OUT;
#构建回溯函数
sub ABC
{
    open OUT,">>","Smith-Waterman.txt";
    my ($i,$j,@data)=@_;
    my @tmp=();
    if (($data[$i-1][$j]-2) == $data[$i][$j])
    {
        push(@tmp,[($i-1),$j]);
        #printf("%d:%d",$i-1,$j);
    }
    if (($data[$i][$j-1]-2) == $data[$i][$j])
    {
        push(@tmp,[$i,($j-1)]);
        #printf("%d:%d",$i,$j-1);
    }
    if (($data[$i-1][$j-1]-3)==$data[$i][$j])
    {
        if ($data[$i][0] ne $data[0][$j])
        {
            push(@tmp,[($i-1),($j-1)]);
            #printf("%d:%d",$i-1,$j-1);
        }
    }
    if (($data[$i-1][$j-1]+9)==$data[$i][$j])
    {
        if ($data[$i][0] eq $data[0][$j])
        {
            push(@tmp,[($i-1),($j-1)]);
            #printf("%d:%d",$i-1,$j-1);
        }
    }
    foreach my $n(0..$#tmp)
    {
        print{OUT}("$tmp[$n][0]:$tmp[$n][1]");
        print{OUT}(" "x10);
    }
    print{OUT}("\n");
    close OUT;
    foreach my $n(0..$#tmp)
    {
        ABC($tmp[$n][0],$tmp[$n][1],@data);
    }
    if ($data[$i][$j]==0)
    {
        return;
    }
}

ABC($row,$col,@matrix);

#输出结果
@array=();
@data=();
push(@array,[1,1]);
push(@array,[2,2]);
open IN,"<","Smith-Waterman.txt";
@data=<IN>;
for ($n=$#data;$n>-1;$n--)
{
    while ($data[$n]=~/(\d+)\:(\d+)/g)
    {
        push(@array,[$1,$2]);
    }
}
foreach $m(0..$#array)
{
    if (($array[$m+1][0]-$array[$m][0])==1 && ($array[$m+1][1]-$array[$m][1])==1)
    {
        print("$DNA_2[$array[$m][0]-1]$DNA_1[$array[$m][1]-1]  ");
    }
    if (($array[$m+1][0]-$array[$m][0])==1 && ($array[$m+1][1]-$array[$m][1])==0)
    {
        print("$DNA_2[$array[$m][0]-1]-  ");
    }
    if (($array[$m+1][0]-$array[$m][0])==0 && ($array[$m+1][1]-$array[$m][1])==1)
    {
        print("-$DNA_1[$array[$m][1]-1]  ");
    }
}
close IN;