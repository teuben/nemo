class Nemo < Formula
  desc "NEMO stellar dynamics toolbox"
  homepage "https://github.com/teuben/nemo"
  url "ftp://ftp.astro.umd.edu/pub/nemo/nemo_4.1.0.tar.gz"
  sha256 "85cc828a96735bdafcf29eb6291ca91bac846579bcef7308536e0c875d6c81d7"

  # depends_on "cmake" => :build

  def install
    system "./configure", "--with-yapp=pgplot",
                          "--prefix=#{prefix}"
    system "make", "build"
  end

  test do
    system "false"
  end
end
