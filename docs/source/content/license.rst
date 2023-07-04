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

License
=======

MIT License

Copyright (c) 2023 Kuan-Hao Chao

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
