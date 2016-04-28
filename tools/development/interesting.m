
inverted = [41.9497, 12.5373, 3.4121];
structure = [41.7910, 11.9025, 2.5392];
initial = [44.0128, 13.6747, 4.2320];


models = [0.05, 0.1, 0.2];

figure
semilogy( models, inverted, 'k*')
hold on
semilogy( models, structure, 'b*')
semilogy( models, initial, 'r*')

legend('inverted', 'only structure', 'source and structure')


figure
plot( models, abs(structure - inverted), 'k*')
hold on
plot( models, abs(structure - initial), 'b*')

legend('structure - inverted', 'structure - initial')