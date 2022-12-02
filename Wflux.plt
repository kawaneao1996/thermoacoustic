#test_3_1.pyの結果をグラフにするgnuplotのスクリプト。コメントアウトは#らしい。

set title "Acoustic field in the resonator, position - Work_flux_{W/m^2}"
set xlabel 'position (m)' ; set ylabel 'work_flux_{W/m^2}'
plot "myfile.txt" using 1:5 with lines title "cal"

#ちなみにresult = 'myfile.txt'の列は[X,Pamp,Uamp,Pphi,W]で並んでいる
