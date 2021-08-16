function [X, time] = smart_timer(f, varargin)
% A timer that records the execution time 'time' of the functoin 'f' of 
% 'varargin' and returns the result X.
  init = tic;
  X = f(varargin{:});
  time = toc(init);
  if time < 1
    timefun = @()(f(varargin{:}));
    time = timeit(timefun);
  end
end