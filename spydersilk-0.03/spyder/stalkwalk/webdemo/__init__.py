header = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <link rel="stylesheet" type="text/css" media="screen, projection" href="spyder.css"/>

  <!-- Foldable menu script -->

  <script type="text/javascript" src="switchcontent.js"></script>
  <script type="text/javascript" src="switchicon.js"></script>
</head>
<body class="brown">
  <div id="col1" class="container">
   <div id="content">

    <div id="maincol">  
"""
footer = """
    <script type="text/javascript">
      var foldmenu=new switchicon("switchgroup1", "div") //Limit scanning of switch contents to just "div" elements
      foldmenu.setHeader('<img src="arrow-down.png"/>','<img src="arrow-up.png"/>') //set icon HTML
      foldmenu.collapsePrevious(false)
      foldmenu.setPersist(true)
      foldmenu.defaultExpanded(0)
      foldmenu.init()
     </script>
    <hr />
   </div>
</body>
<html>
"""

