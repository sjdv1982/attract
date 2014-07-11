$(document).ready(function() {

    // VARIABLES
    var configuration_downloaded = false;
    var uiblock_changed_dict = {};
    var init_clone = ['block-partners-'];

    // FUNCTIONS
    submitForm = function() {

        //Submit the form using jQuery submit. Any pre-submission checking
        //can still be done here.
        $("form").submit();
    }

    trackBlockChange = function() {

        for (var key in uiblock_changed_dict) {
            if (uiblock_changed_dict[key] == 0) {
                $('#nav-' + key).parent().removeClass('changed');
            }
            else {
                $('#nav-' + key).parent().addClass('changed');
            }
        }
    }

    pageIntoView = function($page) {

        var currentID = $page.attr('id');

        //Move corresponding page in view, hide others
        var view = $('#' + currentID.split('-')[1]);
        view.show();
        $('div[id^=page]').each(function() {
            if ($(this)[0] !== view[0]) {
                $(this).hide();
            }
        });

        //Activate menu item, deactivate all other
        $page.addClass('active');
        $('li[id^=view]').each(function() {
            if ($(this)[0] !== $page[0] && $(this).hasClass('active')) {
                $(this).removeClass('active');
            }
        });
    }

    foldAllOthers = function($current) {

        // Fold all similar block to '$current' except current itself
        var currentID = $current.attr('id');
        var rootID = (currentID.split('-').slice(0, 2)).join('-');
        $('div[id^=' + rootID + ']').each(function() {
            if ($(this)[0] !== $current[0] && $(this).hasClass('active')) {
                $(this).find('.level1-container').slideToggle();
                $(this).removeClass('active');
            }
        });
    }

    resetLevelTwo = function($block) {
        $block.find('.level2').each(function() {
            var switch_control = $(this).find('.switch > input');
            if (switch_control.length) {
                if (switch_control.prop('checked')) {
                    $(this).find('.group-container').show();
                }
                else {
                    $(this).find('.group-container').hide();
                }
            }
        });
    }

    reloadFormBlock = function($block) {

        // Reset a part of the form elements to there default state restoring any default values.
        // Applies to all children of the element containing the reset button. Change the form-fields
        // color to the default unchanged version by removing the active class.
        $block.find(':input').each(function() {
            if ($(this).is(':checkbox') || $(this).is(':radio')) {
                $(this).prop('checked', $(this).data('default'));
            }
            else {
                $(this).val($(this).data('default'));
            }
            $(this).parent().removeClass('active');
        });
        resetLevelTwo($block);
    }
    
    getCloneBlocks = function($blockID) {
        
      var regex = new RegExp("^" + $blockID + "[0-9]+$");
      return $('div[id^=' + $blockID + ']').filter(function( index ) {
        return this.id.match(regex);
      });      
    }
    
    renumberFormBlocks = function($root) {
      
        // Renumber form blocks after removal a block to ensure that the remaining
        // form block ID's nicely increment.
        
        $root = $root.substring(1) + '-';
        var counter = 0; //was 1
        getCloneBlocks($root).each(function() {
            
            // DISABLED: If ID equals 0 it is the root, do nothing.
            var curr_id = $(this).attr('id');
            var curr_count = parseInt(curr_id.substring(curr_id.length -1));
            //if (curr_id != $root + 0) { //DISABLED
              
              $(this).attr('id', $root + counter);
              
              $(this).find('[id]').each(function() {
        
                  //Perform the same replace as above
                  var $th = $(this);
                  var newID = $th.attr('id').replace(curr_count, counter);
                  $th.attr('id', newID);
              });
        
              // Find all elements in $clone that have a name, and iterate using each()
              $(this).find('[name]').each(function() {
        
                  //Perform the same replace as above
                  var $th = $(this);
                  var newID = $th.attr('name').replace(curr_count, counter);
                  $th.attr('name', newID);
              });
              
              //Adjust the block title
              var title = $(this).find('.controls .title > h4');
              title.text(title.text().replace(curr_count + 1, counter + 1));
        
              counter += 1;
            //} //
        });
    }
    
    cloneBlock = function($blockID, $max) {

        // Clone a form block with ID '$blockID'. Always take the block with the id '$blockID'-0 as source.
        // After cloning change all the 'id' and 'name' instances in the cloned DOM tree replacing the number
        // '0' with a new number set to the count of all elements matching '$blockID'. Set the cloned block
        // to active and insert it after the last block of the same kind. Activate the blocks 'remove' button
        // After inserting, compare total number of block instances with '$max', if the same, disable the
        // 'add' button.
        count = getCloneBlocks($blockID).length;        

        // Is the user allowed to add more?
        if (count + 1 <= $max) {
            
            // Create your clone, copy including event binders
            var source = $("#" + $blockID + '0');
            var $clone = $(source).clone(true);
            $clone.attr('id', $clone.attr('id').replace('0', count));

            // Find all elements in $clone that have an ID, and iterate using each()
            $clone.find('[id]').each(function() {

                //Perform the same replace as above
                var $th = $(this);
                var newID = $th.attr('id').replace('0', count);
                $th.attr('id', newID);
            });

            // Find all elements in $clone that have a name, and iterate using each()
            $clone.find('[name]').each(function() {

                //Perform the same replace as above
                var $th = $(this);
                var newID = $th.attr('name').replace('0', count);
                $th.attr('name', newID);
            });

            // Remove all highlights (remove embedded resource markings)
            $clone.find('em').remove();
            
            // Activate the blocks remove button
            var $removebutton = $clone.find('#remove-' + $blockID + count);
            if ($removebutton) {
                $removebutton.show();
            }

            // Adjust the block title
            var title = $clone.find('.controls .title > h4');
            title.text(title.text().replace('1', count + 1));

            
            // Reset the blocks values to defaults as they might have been changed in the session
            reloadFormBlock($clone);
            
            // Insert it after the last block with blockID
            $clone.insertAfter(getCloneBlocks($blockID).filter(":last"));
            
            // Set the new block as 'active'
            foldAllOthers($clone);
            $clone.find('.level1-container').slideToggle();
            $clone.addClass('active');

        }

        if (count + 1 == $max) {
            // Set the opacity of the add button to 50% to indicate a disabled button. Button can still fire.
            $('#' + $blockID + 'add').parent().addClass('disabled');
        }
    };

    $.confirm = function(params) {

        // Displays confirmation/alert box as a screen filling overlay. Easy fade in.
        if ($('#confirmOverlay').length) {
            // A confirm is already shown on the page:
            return false;
        }

        var buttonHTML = '';
        $.each(params.buttons,
        function(name, obj) {

            // Generating the markup for the buttons:
            buttonHTML += '<a href="#" class="button ' + obj['class'] + '">' + name + '</a>';

            if (!obj.action) {
                obj.action = function() {};
            }
        });

        var markup = [
        '<div id="confirmOverlay">',
        '<div id="confirmBox">',
        '<h1>', params.title, '</h1>',
        '<p>', params.message, '</p>',
        '<div id="confirmButtons">',
        buttonHTML,
        '</div></div></div>'
        ].join('');

        $(markup).hide().appendTo('body').fadeIn();

        var buttons = $('#confirmBox .button'),
        i = 0;

        $.each(params.buttons,
        function(name, obj) {
            buttons.eq(i++).click(function() {

                // Calling the action attribute when a
                // click occurs, and hiding the confirm.
                obj.action();
                $.confirm.hide();
                return false;
            });
        });
    }

    $.confirm.hide = function() {
        $('#confirmOverlay').fadeOut(function() {
            $(this).remove();
        });
    }


    // CLICK EVENTS
    // Close the application
    $('#close-app').click(function() {
        $.confirm({
            'title': 'Close Application?',
            'message': (configuration_downloaded) ? 'You did not download any configuration yet. <br /> All data will be lost!\
			             Close ATTRACT Online?': 'All data will be lost! Close ATTRACT Online?',
            'buttons': {
                'Yes': {
                    'class': 'blue',
                    'action': function() {
                        self.close();
                    }
                },
                'No': {
                    'class': 'gray',
                }
            }
        });
    });

    // Show or hide the navigation sidepanel
    $('#show-hide-sidepanel').click(function() {

        $(this).toggleClass('button-active');
        var content = $('#content');
        if (content.position().left == 0) {
            content.animate({ left: 185 });
        }
        else {
          content.animate({ left: 0 });
        }

    });

    // Unfold all level 1 blocks of active form category
    $('#unfold-all').click(function() {
        $('li[id^=view-]').each(function() {
            if ($(this).hasClass('active')) {
                var activePageID = $('#' + $(this).attr('id').split('-')[1]);
                activePageID.find('div[id^=block-]').each(function() {
                    if (!$(this).hasClass('active')) {
                        $(this).addClass('active');
                        var container = $(this).find('.level1-container');
                        $(container).slideToggle();
                    }
                });
            }
        });
    });

    // Make all main navigation items clickable. Move corresponding main form view in view, hide others
    $('li[id^=view-]').click(function() {
        pageIntoView($(this));
    });

    // Make all block-* content clickable. If it contains the 'active' class, fold it otherwise unfold.
    $('.controls > .title').click(function() {
        var parent = $(this).closest('div[id^=block-]');
        var container = $(parent).find('.level1-container');
        $(container).slideToggle();
        if ($(parent).hasClass('active')) {
            $(parent).removeClass('active');
        }
        else {
            $(parent).addClass('active');
        }
        foldAllOthers(parent);
    });

    // Show hide level2 blocks on change of the switch
    $('.switch > input').change(function() {
        var parent = $(this).closest('.level2');
        parent.find('.group-container').slideToggle();
    });

    // Check if form fields default value has changed and alter color if so to indicate.
    $(':input').change(function() {
        var blockID = $(this).closest('div[id^=block-]').attr('id');
        if (blockID) {
            var blockID = (blockID.split('-').slice(0, 2)).join('-');
        }

        var changed = false;
        if ($(this).is(':checkbox') || $(this).is(':radio')) {
            if ($(this).prop('checked') !== $(this).data('default')) {
                var changed = true;
            }
        }
        else {
            if ($(this).val() !== $(this).data('default')) {
                var changed = true;
            }
        }

        if (changed) {
            if (blockID in uiblock_changed_dict) {
                uiblock_changed_dict[blockID] += 1
            };
        }
        else {
            if (blockID in uiblock_changed_dict) {
                uiblock_changed_dict[blockID] -= 1
            };
        }
        $(this).parent().toggleClass('active');

        trackBlockChange();
    });

    // Remove a block. Show confirmation alert before removing. Enable 'Add <block name>' button again
    // by removing 'disabled' class if there.
    $('div[id^=remove-block-]').click(function() {

        var blockID = "#" + $(this).attr('id').replace(/^remove-/g, '');
        var rootID = (blockID.split('-').slice(0, -1)).join('-');

        $.confirm({
            'title': 'Delete Confirmation',
            'message': 'You are about to delete this item. <br />It cannot be restored at a later time! Continue?',
            'buttons': {
                'Yes': {
                    'class': 'blue',
                    'action': function() {
                        $(blockID).remove();
                        renumberFormBlocks(rootID);
                        $(rootID + '-add').parent().removeClass('disabled');
                    }
                },
                'No': {
                    'class': 'gray',
                }
            }
        });

    });

    // Reset the full form back to its default value set. do ask for confirmation first
    $('#reload-form').click(function() {
        $.confirm({
            'title': 'Reload Form?',
            'message': 'Are you sure you want to reset all settings back to there default values? this cannot be restored.',
            'buttons': {
                'Yes': {
                    'class': 'blue',
                    'action': function() {
                        reloadFormBlock($('form'));
                        $('div[id^=block-]').each(function() {
                            var blockID = $(this).attr('id');
                            uiblock_changed_dict[(blockID.split('-').slice(0, 2)).join('-')] = 0;
                        });
                        trackBlockChange();
                    }
                },
                'No': {
                    'class': 'gray',
                }
            }
        });
    });

    // Reset a part of the form elements to there default state restoring any default values.
    // Applies to all children of the element containing the reset button. Change the form-fields
    // color to the default unchanged version by removing the active class. Reset the
    // uiblock_changed_dict for the block ID to 0 to indecate that there has been no change made
    // to default values yet.
    $('.reload-form-block').click(function() {
        var curr_block = $(this);
        $.confirm({
            'title': 'Reload this part of the form?',
            'message': 'Are you sure you want to reset this part of the form to its default values? this cannot be restored.',
            'buttons': {
                'Yes': {
                    'class': 'blue',
                    'action': function() {
                        var blockID = curr_block.closest('div[id^=block-]');
                        reloadFormBlock(blockID);
                        uiblock_changed_dict[(blockID.attr('id').split('-').slice(0, 2)).join('-')] = 0;
                        trackBlockChange();
                    }
                },
                'No': {
                    'class': 'gray',
                }
            }
        });
    });


    // TASKS PERFORMED JUST AFTER THE DOCUMENT FINISHED LOADING
    
    // DISABLED: Hide all 'remove' buttons from form blocks with pattern block-*-0 
    //$('div[id^=block-]').each(function() {
    //    if ($(this).attr('id').match('-0$')) {
    //        $('#remove-' + $(this).attr('id')).hide();
    //    }
    //});

    // Store all form defaults in jQuery data objects associated to the form element
    $(':input').each(function() {
        if ($(this).is(':checkbox') || $(this).is(':radio')) {
            $(this).data('default', $(this).prop('checked'));
        }
        else {
            $(this).data('default', $(this).val());
        }
    });

    // Add checked class to all checked checkboxes and radiobuttons.
    // Needed to Firefox to get the active style to work
    $(':input').each(function() {
        if ($(this).is(':checkbox') || $(this).is(':radio')) {
            if ($(this).prop('checked')) {
                $(this).closest('.field-container').addClass('checked');
            }
        }
    });

    // Fold/unfold level 2 form elements
    resetLevelTwo($('form'));

    // Make initial clones for a number of form blocks if they occur in
    // one copy only. Prefilled forms that have multiple copies already are not cloned.
    for (var i in init_clone) {
      if ($('div[id^=' + init_clone[i] + ']').length == 1) {
        cloneBlock(init_clone[i], 10);
      }
    }
  
    // Move the 'Partners' form category into view by default
    pageIntoView($('#view-page1'));

    // Add all UI block names to the uiblock_changed_dict to keep track of changes made.
    $('div[id^=block-]').each(function() {
        var blockID = $(this).attr('id');
        uiblock_changed_dict[(blockID.split('-').slice(0, 2)).join('-')] = 0;
    });

});
