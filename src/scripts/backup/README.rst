CollecTF backup
===============

After a few scary incidents, we finally realized that we must put some
shell scripts that would be run as cron jobs.

This repository contains backup scripts which

- runs ``mysqldump`` to save entire database as a ``.sql`` file every day,
-  copies the latest ``.sql`` file to lab server machine every week.

If something goes wrong, the script sends an email about the issue.

Setup ``ssh`` or ``scp`` login without password
-----------------------------------------------

To do secure copy using ``scp`` without password, local-host’s public
key needs to be saved in remote-host’s ``authorized_keys`` file.

Create public and private keys on local-host
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: shell


    sefa@local-host$ [Note: You are on local-host here]

    sefa@local-host$ ssh-keygen
    Generating public/private rsa key pair.
    Enter file in which to save the key (/home/sefa/.ssh/id_rsa):[Enter key]
    Enter passphrase (empty for no passphrase): [Press enter key]
    Enter same passphrase again: [Press enter key]
    Your identification has been saved in /home/sefa/.ssh/id_rsa.
    Your public key has been saved in /home/sefa/.ssh/id_rsa.pub.
    The key fingerprint is:
    33:b3:fe:af:95:95:18:11:31:d5:de:96:2f:f2:35:f9 sefa@local-host

Copy the public key to remote-host
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

    sefa@local-host$ ssh-copy-id -i ~/.ssh/id_rsa.pub remote-host
    sefa@remote-host's password:
    Now try logging into the machine, with "ssh 'remote-host'", and check in:

    .ssh/authorized_keys

    to make sure we haven't added extra keys that you weren't expecting.

Test: login to remote-host without entering the password
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

    sefa@local-host$ ssh remote-host
    Last login: Sun Nov 16 17:22:33 2008 from 192.168.1.2
    [Note: SSH did not ask for password.]

    sefa@remote-host$ [Note: You are on remote-host here]

Now backup files can be saved to the remote-host using ``scp`` without
password.

Setup cron job
--------------

To run two scripts periodically, they need to be in the crontab file.
Crontab file can be edit by typing

.. code:: example

    crontab -e

which brings the editor to edit the crontab file. For example, to run
one of the scripts every day and the other one every week

::

    0 2 * * * /local/backups/collectf_backup.sh
    0 3 * * 0 /local/backups/scp_to_lab_server.sh
