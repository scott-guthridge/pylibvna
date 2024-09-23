use Graphics::Fig;
use strict;
use warnings;

my $fig = Graphics::Fig->new();

$fig->polygon(12, 3, { rotation => 180 / 12 });
my $bbox = $fig->getbbox();
$fig->translate([-$bbox->[0][0] + 0.5, -$bbox->[0][1] + 0.5]);
$fig->box([[0, 0], [7, 7]]);

$fig->export("dodekagon.png")
