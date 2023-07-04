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

Conclusion
==========

splam is a new deep-learning-based splice junction recognizer that helps us to better understand splice junctions. 


We propose a novel splice junction recognition model called SPLAM, which is a deep convolutional neural
network being able to train on laptop computers. Our approach incorporates several innovative training
concepts. First, we leverage expression information from thousands of samples across different human
tissues in the GTEX project to enhance the quality of splice junctions. Furthermore, we propose a distinctive
method for selecting negative splice junction-like sequences, with a specific focus on those supported by
soley one single splice alignment situated on the opposite strand of the protein coding genes. 