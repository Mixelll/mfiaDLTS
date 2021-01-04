%% Modified by Michael Leitsin for Model 350

function installed = isLakeshoreInstalledTCP()

% Initialize communication to temperature controller.
sess = tcpclient('132.68.54.189',8004);
% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(sess)
    sess = gpib('NI', 0, 12);
else
    fclose(sess);
    sess = sess(1);
end

installed = 0;
r = sess.read;
native2unicode(r)
try
fopen(sess)
fprintf(sess, '*idn?');
pause(.05);

cut = 1:10;
idnCheck = 'LSCI,MODEL350,0,032301';
idn = fscanf(sess);

if strcmp(idn(cut),idnCheck(cut))
    installed = 1;
else
    installed = 0;
end

catch err
    disp('Cannot connect to Lakeshore 350!')
    disp(err.message)
    installed = 0;
end
% Close communication.
fclose(sess)
end