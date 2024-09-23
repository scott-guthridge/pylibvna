#!/usr/bin/perl
use Graphics::Fig;
use Math::Trig;
use strict;
use warnings;

my $N = 3;

my $WIDTH = 7.5;		# overall width of figure
my $DOT_RADIUS = 0.03;

my $DEGREE = 180.0 / pi;
my $ARROW_SPACING = 1.0 / 3.0;

#
# drawDot: draw a filled black dot
#
sub drawDot {
    my $fig   = shift;
    my $point = shift;

    $fig->circle({
	center => $point,
	r => $DOT_RADIUS,
	areaFill => "black",
	depth => 10
    });
}

sub drawArrowLine {
    my $fig       = shift;
    my $points    = shift;
    my $distance  = shift;	# arrow distance from end
    my $label     = shift;
    my $fliplabel = shift;
    my $options   = shift;
    my $offset = 0.062085;	# (inches) put center of arrow on point

    # Uncomment for LaTeX output
    #$label =~ s/([A-Za-z]+)(\d+)/\$\\mathit{$1}_{\\mathrm{$2}}\$/;

    my $dx = $points->[1][0] - $points->[0][0];
    my $dy = $points->[1][1] - $points->[0][1];
    my $length = sqrt($dx*$dx + $dy*$dy);
    my $angle = -atan2($dy, $dx) / pi * 180.0;
    $length -= $distance;

    $fig->begin($options);
    $fig->moveto($points->[0]);
    $fig->lineto($length + $offset, $angle, { arrowMode => "forw" });
    $fig->moveto($points->[0]);
    $fig->moveto($length, $angle);
    $fig->lineto($points->[1]);
    $fig->begin();
    $fig->moveto($points->[0]);
    $fig->moveto($length, 0);
    my $label_off = 0;
    if ((abs($angle) > 90) ^ !!$fliplabel) {
	$label_off = 0.1;
    }
    $fig->moveto(0.085 + $label_off, $fliplabel ? -90 : 90);
    $fig->text($label, { rotation => abs($angle) > 90 ? 180 : 0,
			 textJustification => "center" });
    $fig->moveto($points->[0]);
    $fig->rotate($angle);
    $fig->end();
    $fig->end();
}

#
# drawConnectorArrow: draw an arrow between the error box and DUT
#   $fig:        fig object
#   $start:      reference to [x, y] pair giving start position
#   $startAngle: angle in radians of next point relative to start
#   $end:        reference to [x, y] pair giving end position
#   $endAngle:   angle in radians of previous point relative to end
#   $S:          length of a side of the error box polygon
#   $reverse:    if true, make arrow point from DUT to error box
#   $label:      label for arrow
#
sub drawConnectorArrow {
    my $fig        = shift;
    my $start      = shift;
    my $startAngle = shift;
    my $end        = shift;
    my $endAngle   = shift;
    my $S          = shift;
    my $reverse    = shift;
    my $label      = shift;

    #
    # Find the point in the direction endAngle with distance S from end.
    #
    my $endPrev = [
	$end->[0] + $S * cos($endAngle),
	$end->[1] - $S * sin($endAngle),
    ];

    #
    # Move to start position.
    #
    $fig->begin({ position => $start });
    $fig->moveto($start);

    #
    # Draw outward from the starting point a distance of S, adding
    # the reverse arrow if indicated.
    #
    if ($reverse) {
	$fig->lineto({
	    distance => $ARROW_SPACING * $S,
	    heading  => $startAngle * $DEGREE,
	});
	$fig->lineto({
	   distance  => (1.0 - $ARROW_SPACING) * $S,
	   heading   => $startAngle * $DEGREE,
	   arrowMode => "back"
	});
    } else {
	$fig->lineto({
	    distance => $S,
	    heading  => $startAngle * $DEGREE,
	});
    }

    #
    # Draw to endPrev, first horizontally, then vertically.
    #
    my $cur = $fig->getposition();
    if ($cur->[0] != $endPrev->[0]) {
	$fig->lineto([$endPrev->[0], $cur->[1]]);
    }
    if ($cur->[1] != $endPrev->[1]) {
	$fig->lineto($endPrev);
    }

    #
    # Connect to the DUT.
    #
    if (!$reverse) {
	$fig->lineto({
	   distance  => (1.0 - $ARROW_SPACING) * $S,
	   heading   => $endAngle * $DEGREE + 180.0,
	   arrowMode => "forw"
	});
    }
    $fig->lineto($end);
    $fig->end();
}

#
# Find the number of sides of each polygon.
#
my $calN = 4 * $N;
my $dutN = 2 * $N;

#
# Find the central angles of both polygons in radians.
#
my $calAngle = 2 * pi / $calN;
my $dutAngle = 2 * pi / $dutN;

#
# Find the radius of the error box polygon such that the overall figure
# width is just under $WIDTH.
#
my $calR = $WIDTH * cos(pi / $calN) / 
	(1 + 2 * cos(pi / $dutN) + 4 * sin(pi / $calN));

#
# Find the length of a side of both the error box and DUT polygons.
#
my $S = 2 * $calR * sin(pi / $calN);

#
# Find the radius of the DUT polygon.
#
my $dutR = 0.5 * $calR * sec(pi / $calN);

#
# Find the distance between the center of each polygon and the midpoint
# of a side.
#
my $calB = $calR * cos(pi / $calN);
my $dutB = $dutR * cos(pi / $dutN);

#
# Find the centers of each polygon.
#
my $calX0 = $S + $calB;
my $calY0 = $S + $calB;
my $dutX0 = $calX0 + $calB + 2 * $S + $dutB;
my $dutY0 = $calY0;

#
# Create the drawing object.
#
my $fig = Graphics::Fig->new({
    arrowHeight => "0.1",
    arrowStyle => "filled-triangle",
    arrowWidth => ".1",
    #fontFlags => "+special",		# don't quote LaTeX $S_11$
    lineThickness => "2pt",		# thick lines
    styleVal => ".0375",		# dash spacing
});

#
# Draw the points of each polygon.
#
my @calPoints;
for (my $i = 0; $i < $calN; ++$i) {
    my $angle = $calAngle * ($i + 0.5) + pi / 2.0;
    my $x =  $calR * cos($angle) + $calX0;
    my $y = -$calR * sin($angle) + $calY0;
    push(@calPoints, [$x, $y]);
    &drawDot($fig, [$x, $y]);
}
my @dutPoints;
for (my $i = 0; $i < $dutN; ++$i) {
    my $angle = $dutAngle * ($i + 0.5);
    my $x =  $dutR * cos($angle) + $dutX0;
    my $y = -$dutR * sin($angle) + $dutY0;
    push(@dutPoints, [$x, $y]);
    &drawDot($fig, [$x, $y]);
}

#
# Draw the input and output points of the error box.
#
for (my $i = 0; $i < $N; ++$i) {
    my $angle = pi/2 + (1 + 2*$i) * $calAngle;
    my $dx =  $S * cos($angle),
    my $dy = -$S * sin($angle),
    my $start;
    my $end;

    $end = $calPoints[2*$i];
    $start = [$end->[0] + $dx, $end->[1] + $dy];
    &drawArrowLine($fig, [$start, $end], $S * (1 - $ARROW_SPACING),
	    sprintf("a%d", $i+1), 0, {});

    $start = $calPoints[2*$i + 1];
    $end = [$start->[0] + $dx, $start->[1] + $dy];
    &drawArrowLine($fig, [$start, $end], $S * $ARROW_SPACING,
	    sprintf("b%d", $i+1), 0, {});
}

#
# Draw the directivity terms.
#
for (my $i = 0; $i < $N; ++$i) {
    my $start = $calPoints[2*$i];
    for (my $j = 0; $j < $N; ++$j) {
	my $end = $calPoints[2*$j+1];
	&drawArrowLine($fig, [$start, $end], $S / 2,
		sprintf("Ed%d%d", $j+1, $i+1), 0, { color => "orange" });
    }
}

#
# Draw the transmission terms.
#
for (my $i = 0; $i < $N; ++$i) {
    my $start = $calPoints[2*$i];
    for (my $j = 0; $j < $N; ++$j) {
	my $end = $calPoints[$#calPoints - 2*$j];
	&drawArrowLine($fig, [$start, $end], $S / 2,
		sprintf("Et%d%d", $j+1, $i+1), 0, { color => "green" });
    }
}

#
# Draw the reflection terms.
#
for (my $i = 0; $i < $N; ++$i) {
    my $start = $calPoints[$#calPoints - (2*$i+1)];
    for (my $j = 0; $j < $N; ++$j) {
	my $end = $calPoints[2*$j+1];
	&drawArrowLine($fig, [$start, $end], $S / 2,
		sprintf("Er%d%d", $j+1, $i+1), 0, { color => "blue" });
    }
}

#
# Draw the port mismatch terms.
#
for (my $i = 0; $i < $N; ++$i) {
    my $end = $calPoints[$#calPoints - 2*$i];
    for (my $j = 0; $j < $N; ++$j) {
	my $start = $calPoints[$#calPoints - (2*$j+1)];
	&drawArrowLine($fig, [$start, $end], $S / 2,
		sprintf("Em%d%d", $j+1, $i+1), 0, { color => "purple" });
    }
}

#
# Draw the internal arrows of the DUT.
#
for (my $i = 0; $i < $N; ++$i) {
    my $start = $dutPoints[2*$i];
    for (my $j = 0; $j < $N; ++$j) {
	my $end = $dutPoints[2*$j+1];
	&drawArrowLine($fig, [$start, $end], $S / 2,
		sprintf("S%d%d", $j+1, $i+1), 0, {});
    }
}

#
# Draw the arrows connecting the error box and DUT.
#
for (my $i = 0; $i < $N; ++$i) {
    my $startAngle = pi/2 - (1 + 2*$i) * $calAngle;
    my $endAngle = (1 + 2*$i) * $dutAngle;
    my $start;
    my $end;

    $start = $calPoints[$#calPoints - 2*$i];
    $end = $dutPoints[2*$i];
    &drawConnectorArrow($fig, $start, $startAngle, $end, $endAngle, $S, 0,
    	sprintf("x%d", $i+1));

    $start = $calPoints[$#calPoints - 2*$i - 1];
    $end = $dutPoints[2*$i + 1];
    &drawConnectorArrow($fig, $start, $startAngle, $end, $endAngle, $S, 1,
    	sprintf("y%d", $i+1));
}

$fig->export("hexagon.pdf");
