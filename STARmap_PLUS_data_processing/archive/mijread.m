function image = mijread(filename)
% mijread, version 1.0
%
% im = mijread;
% im = mijread(filename);
% im = mijread(url);
%
% Read the images or stack of images using the ImageJ library. 
% Files must be in TIFF (uncompressed), PNG, GIF, JPEG, DICOM, BMP, PGM or FITS
% format TIFF images can be 1-bit, 8-bits, 16-bits (unsigned), 32-bit
% (real) or RGB color.
% The image could read from the disk (filename) or from a URL (url). If no
% parameter are given to this function, the program open a dialog box to
% choose a image file.

% Example:
% im = mijread;
% im = mijread('http://www.google.com/intl/en_ALL/images/logo.gif');
% 
% Part of the package MIJ
% http://bigwww.epfl.ch/sage/soft/mij/
%
% Installation:
% 1. Put the ij.jar (ImageJ) in the folder java of Matlab
% 2. Put the MIJ.class (MIJ) in the folder java of Matlab
%
opener = ij.io.Opener();
if isjava(opener) == 0
    sprintf('%s', 'MIJ Message: the ImageJ is not properly installed in the java folder of Matlab.')
    image = 0;
    return
end

if (nargin == 0)
    path = pwd;
    dlg = ij.io.OpenDialog('Select an image', '');
    path = dlg.getDirectory();
    if (isjava(path)==0)
        image = 0;
        return
    end
    name = dlg.getFileName();
    if (isjava(name)==0)
        image = 0;
        return
    end
    name = dlg.getDirectory().concat(dlg.getFileName());    
    opener.openImage(name).show();
else
    name = java.lang.String(filename);
    if (name.startsWith('http:') == 1)
         opener.openURL(name).show();
    else
        file = java.io.File(name);
        existence = file.exists();   
        if existence==0;
            sprintf('%s', 'MIJ Message: this file do not exist')
            image = 0;
            return
        end
        opener.openImage(name).show();
    end
    
end

image = MIJ.getCurrentImage();
