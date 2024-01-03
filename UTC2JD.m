%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function JD = UTC2JD(yr, mo, dy, hr, min, sec)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function converts Gregorian date to Julian date through a series of
% parametric calculations.

%--------------------------------------------------------------------------
% Input     - Gregorian Date (YYYY, MM, DD, HH, min, sec)
% Output    - Julian Date (days)
%--------------------------------------------------------------------------

    a = 367 * yr;
    b = fix((7 * (yr + fix((mo + 9) / 12))) / 4);
    c = fix((275 * mo) / 9);
    d = dy;
    e = 1721013.5;
    f = ((sec / 60 + min) / 60 + hr) / 24;
    JD = a - b + c + d + e + f;
end
