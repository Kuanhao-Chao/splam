.. raw:: html

  <script type="text/javascript">
    var observer = new MutationObserver(function(mutations) {
        const dark = document.body.dataset.theme == 'dark';
        console.log(dark);
        document.getElementsByClassName('mainlogo')[0].src = dark ? '../_images/jhu-logo-white.png' : "../_images/jhu-logo-dark.png";
        console.log(document.getElementsByClassName('mainlogo')[0].src);
    })
    observer.observe(document.body, {attributes: true, attributeFilter: ['data-theme']});
    console.log(document.body);
  </script>
  <link rel="preload" href="../_images/jhu-logo-dark.png" as="image">


.. image:: ../image/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, mainlogo
   :align: center

|

Contact
=======

Please use the issue tracker on GitHub for all bug reports. That will help us keep up to date with things.

If you have any new ideas that can improve splam, feel free to contact me <kuanhao.chao@gmail.com>

You can also schedule a `coffee chat with me <https://calendly.com/kuanhao-chao/30min>`_.
