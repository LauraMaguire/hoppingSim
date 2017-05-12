% example model

function [ out ] = model( p, param )

  % assign
  a = p(1);
  b = p(2);
  c = param.c;
  d = param.d;

  % calculate
  out = ( a + b + c + d );

end

