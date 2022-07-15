;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

((c++-mode
  (c-file-style . "linux")
  (c-file-offsets
   (innamespace . nil))
  (c-basic-offset . 4)
  (indent-tabs-mode . nil)
  (fill-column . 100))
 (nil
  (projectile-project-compilation-cmd . "CXXFLAGS=-g make -j4")
  (compile-command . "CXXFLAGS=-g make -j4")))
