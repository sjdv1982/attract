ATTRACT online user interface
=============================


version: 1.0
author : StickyBits (Marc van Dijk), www.stickybits.nl, mvdijk@stickybits.nl
         for Biomolecular Dynamics group, Technische Universitat Munchen
date   : 28-08-2013


Package content
---------------
root:
  - .htaccess: best practice definitions for cross-browser content serving. Needed by ATTRACT online to 
     serve svg files.
  - apple-touch-icon-*: different sized .png images of the ATTRACT icons used by Apple branded mobile
    devices to add as shortcut icon. Unlikely to be used but it is possible.
  - favicon.png: ATTRACT specific browser favicon image
  - index.html: Main application HTML file
  - documentation.html: HTML for ATTRACT application documentation 
  - README.txt: this file :-)
  css:
    - attract.css: ATTRACT online specific Cascade Style Sheet.
    - documentation.css: equals in many respects the attract.css file but adds documentation specific 
      styles and UI form styles are removed.
    - normalize.css: resets major styles to ensure cross-browser compatibility.
    - attract.min.css, documentation.min.css, normalize.min.css: minimized version of the above style
      sheets. Makes them load faster. Made using the online service: http://www.lotterypost.com/css-compress.aspx
  img:
    - All ATTRACT online images all as .svg files. Loads fast and scales resolution independent.
  js:
    - jquery-1.10.2.min.js: compressed version of the jQuery javascript library version 1.10.2
    - modernizr-2.6.2-respond-1.1.0.min.js: browser feature detection library. Not used at the moment.
    - main.js: uncompressed version of ATTRACT UI jQuery functions.
    - main.min.js: compressed version of main.js made using online service: http://jscompress.com
  src:
    - Adobe Illustrator CS6 and SVG source files for ATTRACT UI icons, Logo and initial design mockups

Compatibility
-------------
ATTRACT online UI is tested on:
- FireFox 23.0.x
- Safari 6.0.x
- Internet Explorer 8 and 9
- iPhone and iPad

Known Issues
------------
The ATTRACT online UI is functional in all browsers tested. There are slight differences in visual dynamics and graphical
elements between browsers:
- Internet Explorer does not support CSS3 transitions which makes tooltips and on/of switch element transitions to be 
  instantaneous instead of appearing gradually and with delay (for tooltip)
- input, checkbox and select elements are cannot be fully styled in any of the browsers. Safari (WebKit) provides to most
  support followed by FireFox endng with IE. Specific mozilla definitions in the css allow styling of most elements in 
  FireFox except for the "down arrow" for selection dropdowns. IE adds the inability to style the 'browse' button for files.

HTML5 support
-------------
The interface is made using modern HTML5 definitions. Most of these are not visually apparent but they add to simplicity 
of the code, future ready and (more) cross-browser compatible. In addition, HTML5 offers allot of additional functionality
for form elements such as placeholders for text fields and text areas, client side validation and number fields. The latter
for instance allows to define minimum and maximum values and step-size with an up/down click button to change the value.
Number fields currently work in Safari but not in FireFox and IE where just act as regular text fields. I nevertheless added
the functionality to the interface.Same is true for placeholders. I did not add client side validation but it might be 
possible although for simple cases only.

UI features
-----------
- Application like interface design, scales to fill the browser window. Header and sidebar are fixed, overflow is made
  scrollable
- The header menu bar allows for (from left to right): show/hide the left side form catagory menu; reset all form elements
  to there default values; unfold all foldable UI blocks in the current view; a small dropdown menu with custom links to
  contact information, documentation e.d.; main application quite button
- UI actions that result in changes that cannot be undone (reset, delete, quite) show a confirmation dialog first
- Main form catagories are shown in the left sidebar menu. Active ones are colored. They each show a div element containing
  corresponding form elements while hiding others (jQuery controlled). These div's are of class 'category-container' and id
  'page1' to 'page9'.
- When any of the form elements are changed to a state with a value other than the default the element in colored blue. A
  blue dot next to the form catagory item in the left side menu indicates that values have been set/changed for that category.
- There are reset buttons for each foldable form block which enable the user to reset the form elements for that block only.
- Top level foldable form blocks (level1) can have several second level foldable blocks (level2) that are controlled using
  an on/off switch. Often these switched set checkbox values that state rather on not to use a particular ATTRACT feature
  displaying the corresponding fields ones switched to the 'on' state.
- Form elements that allow for multiple copies like multiple docking partners only have to present onces in the initial form.
  This root (or source) block with id 'block-<catagory name>-number' is cloned on the client side when, for example, multiple
  partners are requested. The 'add partner' button performs this action. It takes as a argument the id of the source block and the
  maximum allowed number of copies. When the maximum is reached the button is disabled. It's enabled again when the user 
  removes a block. Specificly for partners the source element is cloned once on application initiation to always start with
  two partners at least.
- The 'Get configuration' button is the form submit function. It's controlled by the javascript equivalent. Because it is a
  function you might do some pre-submission checks if you like to.
- Hovering over a form element for 1.5 sec shows a tooltip which contains a link to corresponding documentation. The link points
  to documentation.html and then navigates to the section describing the feature using document anchors.
  
Code notes
----------
- All default HTML form elements are wrapped in a number of div's to enable styling, dynamic behaviour and tooltips. Make
  sure to always use this setup otherwise it breaks the interface. If you are unsure on how the HTML functions please contact
  me.

HTML structure
--------------
Description of the HTML hierarchy with respect to the different form levels, toplevel down. All form levels are contained within the
<form></form> element.  

-  Div of class "category-container": top level wrapper around all category form elements like "Partners", "Grids" etc. Each of 
   these divs has a unique id ad "page<number>". This id is used to enable jQuery driven switching between category (or pagers) by
   clicking on the category name in the left side menu.
   
   - Div of class "form-category-description": small block containing a short description of the category and optionally a block (div of
     class "category-controls") with buttons that allow for duplication (cloning) of particular form blocks within the category       
   - Div of class "form-category": simple wrapper div containing all form elements for the given category
    
      - Div of class "level1": first level collection of form elements like "Partner 1". Level 1 is a foldable block with header containing
        a title, folding icon, reload button and remove button (all contained in a div of class "controls"). The form elements themselfs 
        are containg within a wrapper of class "level1-container".
        
        - Div of class "level2": subcategory of form elements with a header title and an optional "slide switch" to fold/unfold the form
          elements belonging to the level2.
            
            - Div of class "group-container": wrapper containing form elements of level2. It's function is pure positioning.
            - Div of class "group-inline": wrapper to enable "inline" positioning of field-containers.
            - Div of class "field-container": this is the main container div around the basic input and select form elements.
            
              Field container:
              - Div of class "title": optional title of the form element. positioned inline to the right of the input element. A field-container
                is allowed to only contain a title for which the text may be wrapped in <p> or <h#>. This is for instance done with
                the headers if Level 2 blocks. 
              - Div of class "tooltip": optional tooltip for the form element. Contains a short description optionally linking to 
                a documentation page. tooltips appear if hovered over input element for 2 seconds.
              - Div of class "field-item": wrapper around the input or select element enabling ATTRACT specific styling. Differentiation in
                styles are enabled by adding additional classes to field-item like "text" "select" "file" "number" or "checkbox". 
              - NOTE: if the field-container has class "switch" added, it result in the blue sliding switcher displayed that will enable 
                folding/unfolding of level2 blocks. It is required that the field container containing the switch field-container is 
                a child of level2. The switch field-container may contain a input element of type checkbox to couple switching to a
                persistent choice. It's not required though to enable switch based folding/unfolding.