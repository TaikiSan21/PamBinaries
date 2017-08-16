function len = getClickDataLen(click)
% get the length of this click when it's written to file
len = 8+4+8+4+4+2+2+4+4*length(click.delays)+2+4*length(click.angles)+2+4*length(click.angleerrors)+2+4+click.duration;