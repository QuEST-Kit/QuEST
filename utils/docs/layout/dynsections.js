/*
 Edit by Tyson Jones, 10th June 2021
 I don't know why I'm adding this file, but Doxygen doc gen on Ubuntu 18
 is producing html which refers to this file but not the file itself. 
 I've copied it from a MacOS build which does build it, but doesn't seem to 
 render the HTML properly. I've removed reference to a non-existent (at least
 in my Ubuntu build) powerTip() function.
 Finally, I added a hack in $(document.ready) to delete HTML elements which 
 aren't rendering properly, despite my best 4am efforts. 
 Such is my wrath ¯\_(ツ)_/¯
 
 
 @licstart  The following is the entire license notice for the
 JavaScript code in this file.

 Copyright (C) 1997-2017 by Dimitri van Heesch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

 @licend  The above is the entire license notice
 for the JavaScript code in this file
 */
function toggleVisibility(linkObj)
{
 var base = $(linkObj).attr('id');
 var summary = $('#'+base+'-summary');
 var content = $('#'+base+'-content');
 var trigger = $('#'+base+'-trigger');
 var src=$(trigger).attr('src');
 if (content.is(':visible')===true) {
   content.hide();
   summary.show();
   $(linkObj).addClass('closed').removeClass('opened');
   $(trigger).attr('src',src.substring(0,src.length-8)+'closed.png');
 } else {
   content.show();
   summary.hide();
   $(linkObj).removeClass('closed').addClass('opened');
   $(trigger).attr('src',src.substring(0,src.length-10)+'open.png');
 }
 return false;
}

function updateStripes()
{
  $('table.directory tr').
       removeClass('even').filter(':visible:even').addClass('even');
}

function toggleLevel(level)
{
  $('table.directory tr').each(function() {
    var l = this.id.split('_').length-1;
    var i = $('#img'+this.id.substring(3));
    var a = $('#arr'+this.id.substring(3));
    if (l<level+1) {
      i.removeClass('iconfopen iconfclosed').addClass('iconfopen');
      a.html('&#9660;');
      $(this).show();
    } else if (l==level+1) {
      i.removeClass('iconfclosed iconfopen').addClass('iconfclosed');
      a.html('&#9658;');
      $(this).show();
    } else {
      $(this).hide();
    }
  });
  updateStripes();
}

function toggleFolder(id)
{
  // the clicked row
  var currentRow = $('#row_'+id);

  // all rows after the clicked row
  var rows = currentRow.nextAll("tr");

  var re = new RegExp('^row_'+id+'\\d+_$', "i"); //only one sub

  // only match elements AFTER this one (can't hide elements before)
  var childRows = rows.filter(function() { return this.id.match(re); });

  // first row is visible we are HIDING
  if (childRows.filter(':first').is(':visible')===true) {
    // replace down arrow by right arrow for current row
    var currentRowSpans = currentRow.find("span");
    currentRowSpans.filter(".iconfopen").removeClass("iconfopen").addClass("iconfclosed");
    currentRowSpans.filter(".arrow").html('&#9658;');
    rows.filter("[id^=row_"+id+"]").hide(); // hide all children
  } else { // we are SHOWING
    // replace right arrow by down arrow for current row
    var currentRowSpans = currentRow.find("span");
    currentRowSpans.filter(".iconfclosed").removeClass("iconfclosed").addClass("iconfopen");
    currentRowSpans.filter(".arrow").html('&#9660;');
    // replace down arrows by right arrows for child rows
    var childRowsSpans = childRows.find("span");
    childRowsSpans.filter(".iconfopen").removeClass("iconfopen").addClass("iconfclosed");
    childRowsSpans.filter(".arrow").html('&#9658;');
    childRows.show(); //show all children
  }
  updateStripes();
}


function toggleInherit(id)
{
  var rows = $('tr.inherit.'+id);
  var img = $('tr.inherit_header.'+id+' img');
  var src = $(img).attr('src');
  if (rows.filter(':first').is(':visible')===true) {
    rows.css('display','none');
    $(img).attr('src',src.substring(0,src.length-8)+'closed.png');
  } else {
    rows.css('display','table-row'); // using show() causes jump in firefox
    $(img).attr('src',src.substring(0,src.length-10)+'open.png');
  }
}
/* @license-end */

$(document).ready(function() {
  $('.code,.codeRef').each(function() {
    $(this).data('powertip',$('#a'+$(this).attr('href').replace(/.*\//,'').replace(/[^a-z_A-Z0-9]/g,'_')).html());
    
    // @TysonJones
    // removing this reference to undefined function powerTip (but I don't know why)
    // $(this).powerTip({ placement: 's', smartPlacement: true, mouseOnToPopup: true });
  });
  
  // @TysonJones
  // deleting all memtitle elements (appearing above each function doc), which weren't 
  // rendering properly in bootstrap & smartmenu. To future potential employers who
  // discover my crimes here: eat my shorts
  var glitchytitles = document.getElementsByClassName("memtitle");
  while(glitchytitles.length > 0) {
    glitchytitles[0].parentNode.removeChild(glitchytitles[0]);
  }
  
  // @TysonJones
  // Here, we can increase the font size of all function signatures, to make up for the above 
  // deleted memtitles. However, I presently prefer the signature size being the same 
  // as the doc main text
  /*
  var signatureElems = document.getElementsByClassName("memname");
  for (const elem of signatureElems) {
    // main doc text size appears to be 14px
    elem.style.fontSize = "16px";
  }
  */
});
