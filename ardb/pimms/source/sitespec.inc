* This is the Site Specific include file used by PIMMS
* Please change the directory specification below to suit the
* implementation at your site.

        character*256 ddir_name, mdir_name

* VAX/VMS examples
*        data ddir_name / '$disk1:[mukai.nra.pimms.develop.data]' /
*        data mdir_name / '$disk1:[mukai.nra.pimms.develop.models]' /

* UNIX examples
        data ddir_name /
     & '/data/fubar/SCAR/ardb/pimms/data' /
        data mdir_name /
     & '/data/fubar/SCAR/ardb/pimms/models' /
