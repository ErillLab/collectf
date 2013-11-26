(require 'org-publish)
(setq org-publish-project-alist
      '(("collectf"
         :base-directory "."
         :base-extension "org"
         :publishing-directory "build"
         :recursive t
         :html-postamble nil
         :section-numbers nil
         :publishing-function org-publish-org-to-html
         :style "<link rel=\"stylesheet\" href=\"worg.css\" type=\"text/css\" />"
         )))
