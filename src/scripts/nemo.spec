Summary: A Stellar Dynamics Toolbox
Name: nemo
Version: 3.2.2
Release: 1
License: GPL
Group: Development/Tools
URL: http://www.astro.umd.edu/nemo
Source0: ftp://ftp.astro.umd.edu/pub/nemo/%{name}-%{version}.tar.gz

BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildArch: noarch

%description
NEMO is an extendible  Stellar Dynamics Toolbox, following the
Open-Source Software model. It has various programs to create,
integrate, analyze and visualize N-body and SPH like systems, following
the pipe and filter architecture. In addition there are various tools to
operate on images, tables and orbits, including FITS files to
export/import to/from other astronomical data reduction packages. A
large growing fraction of NEMO has been contributed by a growing list of
authors. The source code consist of a little under 1000 files and
150,000 lines of code, mostly C, and some C++ and Fortran. 

%prep

%build

%install
rm -rf $RPM_BUILD_ROOT

mkdir -p $RPM_BUILD_ROOT/etc
touch $RPM_BUILD_ROOT/etc/empty-file

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,-)
/etc/empty-file

