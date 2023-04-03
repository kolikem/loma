#!/bin/sh
chmod +x loma
mv -i loma /usr/local/bin
chmod +x loma_src/*.py
mv -i loma_src /usr/local/bin
echo "Setup is done."
echo "User can use the command 'loma' now."
