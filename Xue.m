% importData;
NameMatch=zeros(31290,1);
for i=1:31290
if ~isempty(findstr(DisplayName{i}, 'xiaojingzi_xue'))
NameMatch(i)=1;
end
end
plot(NameMatch);
for i=1:numel(Date)
Date{i}=Date{i}(1:10);
end
DV  = datevec(Date);  % [N x 6] array
DV  = DV(:, 1:3);   % [N x 3] array, no time
DV2 = DV;
DV2(:, 2:3) = 0;    % [N x 3], day before 01.Jan
Result = cat(2, DV(:, 1), datenum(DV) - datenum(DV2));
figure
histogram(Result(:,2),'BinMethod','integers')
Message=Message(logical(NameMatch));
