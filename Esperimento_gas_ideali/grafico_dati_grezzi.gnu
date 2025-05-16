set grid
set key inside left top box 
set title "Grafico V-1/P compressione" font "Arial,13"
set xrange [0.2:1.2]
set yrange [5:25]
set xlabel "1/P [cm²/kg]" font "Arial,13"
set ylabel "V [cm³]" font "Arial,13"
set tics font ",12"  # Applica la dimensione del font a entrambi gli assi
set xtics 0.01      # Tacche ogni 2 unità sull'asse X
set ytics 0.1  
plot "0gradicompressione.txt" using (1/$1):2 with lines linecolor "red" title "0°C", "15gradicompressione.txt" using (1/$1):2 with lines linecolor "orange" title "15°C", "25gradicompressione.txt" using (1/$1):2 with lines linecolor "pink" title "25°C", "35gradicompressione.txt" using (1/$1):2 with lines linecolor "green" title "35°C", "45gradicompressione.txt" using (1/$1):2 with lines linecolor "blue" title "45°C","55gradicompressione.txt" using (1/$1):2 with lines linecolor "purple" title "55°C", "0gradidilatazione.txt" using (1/$1):2 with lines linecolor "red" title "0°C", "15gradidilatazione.txt" using (1/$1):2 with lines linecolor "orange" title "15°C", "25gradidilatazione.txt" using (1/$1):2 with lines linecolor "pink" title "25°C", "35gradidilatazione.txt" using (1/$1):2 with lines linecolor "green" title "35°C", "45gradidilatazione.txt" using (1/$1):2 with lines linecolor "blue" title "45°C","55gradidilatazione.txt" using (1/$1):2 with lines linecolor "purple" title "55°C"
