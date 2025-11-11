function results = gifti_runtests(pth)

if ~nargin, pth = fileparts(mfilename('fullpath')); end

d = dir(pth);
d([d.isdir]) = [];
d(arrayfun(@(x) isempty(regexp(x.name,'^test_.*\.m$','once')),d)) = [];

results = struct('Passed',{},'Failed',{},'Incomplete',{},'Duration',{});
for i=1:numel(d)
    results(i).Failed = false;
    results(i).Passed = false;
    results(i).Incomplete = false;
    tstart = tic;
    try
        fprintf('%s',d(i).name(1:end-2));
        feval(d(i).name(1:end-2));
        results(i).Passed = true;
    catch
        results(i).Failed = true;
        err = lasterror;
        fprintf('\n%s',err.message);
    end
    results(i).Duration = toc(tstart);
    fprintf('\n');
end

if ~nargout
    fprintf(['Totals (%d tests):\n\t%d Passed, %d Failed, %d Incomplete.\n' ...
        '\t%f seconds testing time.\n\n'],numel(results),nnz([results.Passed]),...
        nnz([results.Failed]),nnz([results.Incomplete]),sum([results.Duration]));
end
