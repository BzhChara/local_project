#````1````
$dna=('CGCTCGGGGGGCATGTCAAATGAGCTCAACAACATCATCTCTCGGACCACAGATGGCGTCTATGAGGGCGTGGCCATCGGCGGGGACAGGTACCCTGGCTCAACATTTATGGATCACGTGTTACGTTACCAGGATACTCCAGGAGTAAAAATGATTGTGGTTCTTGGAGAGATAGGGGGCACTGAGGAATACAAGATCTGCCGGGGCATCCAGGAGGGCCGCCTCACCAAGCCCGTGGTCTGCTGGTGCATCGGGACATGTGCCACCATGTTCTCTTCTGAGGTACAGTTCGGCCATGCCGGAGCTTGCGCCAACCAGGCTTCCGAAACTGCAGTAGCCAAGAACCAGGCTTTGAAGGAGGCAGGAGTGTTTGTGCCCCCGAGCTTTGATGAACTTGGAGAAATCATCCAGTCTGTGTATGAAGATCTTGTGGCCAGAGGAGTCATTGTCCCTGCTCAGGAGGTGCCGCCTCCAACCGTGCCCATGGACTACTCCTGGGCCAGGGAGCTGGGTTTGATCCGCAAACCTGCCTCATTCATGACCAGCATCTGTGACGAGCGAGGACAGGAGCTCATCTATGCGGGCATGCCCATCACCGAGGTCTTCAAGGAGGAGATGGGCATTGGTGGGGTCCTTGGCCTCCTGTGGTTCCAGAGAAGGTTGCCCAAGTATGCCTGCCAGTTCATTGAGATGTGCCTGATGGTGACGGCAGATCACGGGCCAGCTGTGTCTGGGGCTCACAACACCATCATCTGCGCTCGGGCTGGGAAGGACCTGGTTTCCAGCCTCACCTCGGGGCTGCTCACTATTGGGGACCGGTTTGGGGGTGCCCTGGATGCTGCTGCCAAGATGTTCAGCAAAGCCTTTGACAGTGGTATTATCCCCATGGAGTTTGTGAACAAGATGAAGAAGGAAGGAAAGCTG');
@type=qw/AA AG AC AT GA GG GC GT CA CG CC CT TA TG TC TT/;
@DNA=split('',$dna);
$length_1=@type;
$length_2=@DNA;
@num=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
foreach $n(0..$length_1-1)
{
    foreach $m(0..$length_2-1)
    {
        if (substr($dna,$m,2) eq $type[$n])
        {
            $num[$n]++;
        }
    }
}
$total=0;
foreach $n(0..$length_1-1)
{
    $total=$total+$num[$n];
}
foreach $n(0..$length_1-1)
{
    printf("$type[$n]=$num[$n],percent=%f%\n",$num[$n]/$total*100);
}
print("\n");

#````2````
@array=([8,5,0,30,0,0],[0,0,5,0,71,0],[2,32,1,1,0,0],[55,45,42,37,18,15]);
@sum=(0,0,0,0);
foreach $n(0..3)
{
    foreach $m(0..5)
    {
        $sum[$n]=$sum[$n]+$array[$n][$m];
    }
}
print("@sum\n");
@new=();
foreach $i(0..3)
{
    push(@new,@{$array[$i]});
}
@new=sort{$b<=>$a}@new;
print("@new\n");
print("\n");

#````3````
@t=([0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]);
foreach $n(0..3)
{
    foreach $m(0..5)
    {
        $t[$m][$n]=$array[$n][$m];
    }
}
foreach $i(0..$#t)
{
    print("@{$t[$i]}\n");
}
print("\n");

#````4`````
$data=('A mixture of novel glycopeptides from
glycosylation between cold water fish skin gelatin 
hydrolysates and glucosamine (GlcN) via transglutaminase (TGase), 
as well as glycation between fish gelatin hydrolysate and GlcN were 
identified by their pattern of molecular distribution using MALDI-TOF-MS.');
@DATA=split('\n',$data);
foreach $i(0..$#DATA)
{
    @word=split(' ',$DATA[$i]);
    print("$word[$i]\n");
}
print("\n");

#````5````
@DATA=split(' ',$data);
$word=join('*',@DATA);
print("$word\n");
print("\n");

#````6````
@num=(63,96,70,0,9,50,100);
for ($j=0;$j<scalar(@num)-1;$j++)
{
    for ($i=0;$i<scalar(@num)-1-$j;$i++)
    {
        if ($num[$i]>$num[$i+1])
        {
            ($num[$i],$num[$i+1])=($num[$i+1],$num[$i]);
        }
    }
}
print("@num\n");
print("\n");

#````7````
$array=('                      x_cor    y_cor   z_cor
ATOM      1  N   VAL A   1     101.601  38.534  -1.962  1.00 53.29           N  
ATOM      2  CA  VAL A   1     103.062  38.513  -2.159  1.00 47.99           C ');
@data=split('\n',$array);
@N=split(' ',$data[1]);
@C=split(' ',$data[2]);
$distance=(($N[6]-$C[6])**2+($N[7]-$C[7])**2+($N[8]-$C[8])**2)**(1/2);
print("$distance\n");
print("\n");

#````8````
$old=('aagcttgctt tcattagaaa gacgagacag cagctttcca aagatacaca cagcacttga');
$new=('aagcttgctt tcatttgaaa gacgagacag cagctttcca aagattcaca cagcacatga');
@old=split(' ',$old);
$old=join('',@old);
@new=split(' ',$new);
$new=join('',@new);
foreach $i(0..length($old)-1)
{
    if (substr($old,$i,1) ne substr($new,$i,1))
    {
        $n=substr($old,$i,1);
        $m=substr($new,$i,1);
        printf("%d,old:$n,new:$m\n",$i+1);
    }
}
print("\n");

#````9````
@array=(1,3,6,3,7,34,6,3,7,9,3,0,1,4,8,5,55,90,11,24,67,90,39,2,90,34);
@new=();
$n=0;
push(@new,$array[0]);
for ($i=1;$i<scalar(@array);$i++)
{
    for ($j=$n=0;$j<scalar(@new);$j++)
    {
        if ($new[$j] != $array[$i])
        {
            $n++;
        }
    }
    if ($n == scalar(@new))
    {
        push(@new,$array[$i]);
    }
}
print("@new\n");