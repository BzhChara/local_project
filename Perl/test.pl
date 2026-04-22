use strict;
use warnings;

# 转移概率矩阵
my @A = (
    [0.99, 0.02 / 3, 0.01 / 3],
    [0.01 / 2, 0.99, 0.01 / 2],
    [0.02 / 3, 0.01 / 3, 0.99]
);

# 输出概率矩阵
my @B = (
    [0.03, 0.97],
    [0.5, 0.5],
    [0.97, 0.03]
);

# 初始状态分布
my @initp = (0.4, 0.2, 0.4);

# 状态数和状态名
my $S = 3;
my @S_name = ('1', '2', '0');

# 观察序列
my $Y = '1111111011111101111111000000001000000001000001000000000000001011001010101010111111101111111111101111111011111';

# viterbi实现
# 向前计算
my @Probability;
my @Position;
push @Position, [0, 0, 0];

my @first = map { $initp[$_] * $B[$_][substr($Y, 0, 1)] } 0..$S-1;
push @Probability, \@first;

for my $i (1..length($Y)-1) {
    my @ProbabilityTemp;  # 用来存放每一层观察值对应的三种状态概率值
    my @GetPosition;
    for my $j (0..$S-1) {
        my @temp = map { $Probability[$i-1][$_] * $A[$j][$_] } 0..$S-1;
        push @ProbabilityTemp, (max(@temp)) * $B[$j][substr($Y, $i, 1)];
        push @GetPosition, (max_idx(@temp));
    }
    push @Probability, \@ProbabilityTemp;
    push @Position, \@GetPosition;
}

my $traceback_idx = max_idx(@{$Probability[-1]});
my @path = ($S_name[$traceback_idx]);

for my $i (reverse 0..length($Y)-2) {
    push @path, $S_name[$traceback_idx];
    $traceback_idx = $Position[$i][$traceback_idx];
}

# 输出结果
print(join('', reverse @path), "\n");

sub max {
    my $max = shift;
    for (@_) {
        $max = $_ if $_ > $max;
    }
    return $max;
}

sub max_idx {
    my $max_index = 0;
    for my $i (1..$#_) {
        $max_index = $i if $_[$i] > $_[$max_index];
    }
    return $max_index;
}
