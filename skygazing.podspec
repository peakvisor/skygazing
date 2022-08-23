Pod::Spec.new do |s|
  s.name         = "skygazing"
  s.version      = "1.0.0"
  s.summary      = "Lightweight header-only C++ library for astronomical calculations."
  s.homepage     = "https://github.com/peakvisor/skygazing"
  s.license      = { :type => "GPL-3.0", :file => "LICENSE" }
  s.author             = { "eofom" => "eeofom@gmail.com" }
  s.source       = { :git =>  s.homepage , :commit => "62543fc3f5"}
  s.source_files = "lib/*.h", "include/*.h"
end