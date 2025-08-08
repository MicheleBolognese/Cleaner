within Examples;
block prova
    extends .Modelica.Blocks.Icons.Block;
constant Real [2,3] a(start = [1.0,2.0,3.0;5.0,6.0,7.0]);
constant Real [5] b(start = {1.0,1.0,1.0,1.0,1.0});
Real result2[5,2];
parameter Real [4] L={2,6,8,10};

final parameter Real [4] z={sum(L[1:i]) + 0.5*L[i + 1] for i in 0:4 - 1};


equation
for i in 1:5 loop
result2[i,:] = sum(a[:,j] for j in 1:3)*b[i];

end for;

end prova;
