#test_3_1.pyの結果をグラフにするgnuplotのスクリプト。コメントアウトは#らしい。

set title "Acoustic field in the resonator, position - P_{amp}"
set xlabel 'position (m)' ; set ylabel 'P_{amp}'
plot "myfile.txt" using 1:2 with lines title "cal"

#ちなみにresult = 'myfile.txt'の列は[X,Pamp,Uamp,Pphi,W]で並んでいる
