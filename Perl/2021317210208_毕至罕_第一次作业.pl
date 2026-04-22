#1
printf ("*****1*****\n");
chomp ($name=<STDIN>);
printf ("$name, welcome to Perl programming course!\n");
printf ("\n");
#2
printf ("*****2*****\n");
printf ("Enter the number:");
chomp ($input_1=<STDIN>);
printf ("%d,%X,%o\n",$input_1,$input_1,$input_1);
printf ("\n");
#3
printf ("*****3*****\n");
printf ("Enter two numbers:");
chomp ($input_2=<STDIN>);
chomp ($input_3=<STDIN>);
$data_1=$input_2+$input_3;
$data_2=$input_2-$input_3;
$data_3=$input_2*$input_3;
$data_4=$input_2/$input_3;
$data_5=$input_2**$input_3;
printf ("$input_2+$input_3=$data_1
$input_2-$input_3=$data_2
$input_2*$input_3=$data_3
$input_2/$input_3=$data_4
$input_2**$input_3=$data_5\n");
printf ("\n");
#4
printf ("*****4*****\n");
$PI=3.141592654;
printf ("Enter the number:");
chomp ($input_4=<STDIN>);
$C=2*$PI*$input_4;
$S=$PI*($input_4**2);
printf ("$C,$S\n");
printf ("\n");
#5
printf ("*****5*****\n");
printf ("Enter the number:");
chomp ($input_5=<STDIN>);
if ($input_5>50)
{
    printf ("The data is bigger than 50!\n");
}
else
{
    printf ("The data is less than 50!\n");
}
printf ("\n");
#6
printf ("*****6*****\n");
$n=1;
$num=1;
$sum=0;
while ($num<=50)
{
    printf ("$num ");
    $sum=$sum+$num;
    $num=$num+2;
    $n=$n+1;
    if ($n>5)
    {
        $n=1;
        printf ("\n");
    }
}
printf ("sum is:$sum\n");
printf ("\n");
#7
printf ("*****7*****\n");
$str_1='abcdfeg';
$str_2='higklmn';
while ($str_1)
{
    $str=chop ($str_1);
    printf ("$str");
}
printf ("\n");
while ($str_2)
{
    $str=chop ($str_2);
    printf ("$str");
}
printf ("\n\n");
#8
printf ("*****8*****\n");
printf ("Please input your ID:");
$n=0;
$name='perlbaby';
$passwaord='helloperl';
chomp ($input_6=<STDIN>);
if ($input_6 eq $name)
{
    printf ("Welcome perlbaby!\n");
    printf ("Please input your password:");
    while (1)
    {
        chomp ($input_7=<STDIN>);
        if ($input_7 eq $passwaord)
        {
            printf ("Your password is correct!\n");
            last;
        }
        else
        {
            printf ("Your password is wrong!\n");
            printf ("Please input again:");
            chomp ($input_7=<STDIN>);
            $n=$n+1;
            if ($n==3)
            {
                printf ("Sorry, please try again tomorrow!\n");
                last;
            }
        }
    }
}
else
{
    printf ("Your input ID is not authorized!\n");
}
printf ("\n");
#9
printf ("*****9*****\n");
$n=$m=0;
@str=qw/* * * * * */;
while ($n<=5)
{
    while ($m<=$n)
    {
        printf ("$str[$m]");
        $m=$m+1;
    }
    printf ("\n");
    $n=$n+1;
    $m=0;
    if ($n==3)
    {
        $n=$n+1;
        next;
    }
}
printf ("\n");
#10
printf ("*****10*****\n");
$protein='MNAPERQPQPDGGDAPGHEPGGSPQDELDFSILFDYEYLNPNEEEPNAHKVASPPSGPAYPDDVLDYGLKPYSPLASLSGEPPGRFGEPDRVGPQKFLSAAKPAGASGLSPRIEITPSHELIQAVGPLRMRDAGLLVEQPPLAGVAASPRFTLPVPGFEGYREPLCLSPASSGSSASFISDTFSPYTSPCVSPNNGGPDDLCPQFQNIPAHYSPRTSPIMSPRTSLAEDSCLGRHSPVPRPASRSSSPGAKRRHSCAEALVALPPGASPQRSRSPSPQPSSHVAPQDHGSPAGYPPVAGSAVIMDALNSLATDSPCGIPPKMWKTSP';
$length=length($protein);
printf ("length:$length\n");
$n=0;
while ($protein)
{
    $str=chop ($protein);
    if ($str eq 'A')
    {
        $n=$n+1;
    }
}
printf ("proportion:%d%",$n/$length*100);
printf ("\n");