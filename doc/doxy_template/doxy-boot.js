// LICENSE: This file (but not the other files in this directory) is
// licensed by the Apache License Version 2.0.
// See the file doxy-boot.js.LICENSE for details.
//
// This file is taken from
// https://github.com/Velron/doxygen-bootstrapped
// and some modification and extensions are applied!

$( document ).ready(function() {

  $("ul.tablist").addClass("nav nav-pills nav-justified");
  $("li.current").addClass("active");

  $("#nav-path > ul").addClass("breadcrumb");

  $("table.params").addClass("table");
  $("div.ingroups").wrapInner("<small></small>");
  $("div.levels").css("margin", "0.5em");
  $("div.levels > span").addClass("btn btn-default btn-xs");
  $("div.levels > span").css("margin-right", "0.25em");

  $("table.directory").addClass("table table-striped table-hover table-bordered table-condensed");
  $("div.summary > a").addClass("btn btn-default btn-xs");
  $("table.fieldtable").addClass("table");
  $(".fragment").addClass("well");
  $(".memitem").addClass("panel panel-default");
  $(".memproto").addClass("panel-heading");
  $(".memdoc").addClass("panel-body");
  $("span.mlabel").addClass("label label-success");

  $("table.memberdecls").addClass("table");
  $("[class^=memitem]").addClass("active");

  $("div.ah").addClass("btn btn-default");
  $("span.mlabels").addClass("pull-right");
  $("table.mlabels").css("width", "100%")
  $("td.mlabels-right").addClass("pull-right");

  $("div.ttc").addClass("panel panel-primary");
  $("div.ttname").addClass("panel-heading");
  $("div.ttname a").css("color", 'white');
  $("div.ttdef,div.ttdoc,div.ttdeci").addClass("panel-body");

  $('#MSearchBox').parent().remove();

  $('div.fragment.well div.line:first').css('margin-top', '15px');
  $('div.fragment.well div.line:last').css('margin-bottom', '15px');
	
  $('table.doxtable').removeClass('doxtable').addClass('table table-striped table-bordered').each(function() {
    $(this).prepend('<thead></thead>');
    $(this).find('tbody > tr:first').prependTo($(this).find('thead'));
  
    $(this).find('td > span.success').parent().addClass('success');
    $(this).find('td > span.warning').parent().addClass('warning');
    $(this).find('td > span.danger').parent().addClass('danger');
  });
	
	

  if($('div.fragment.well div.ttc').length > 0)
  {
    $('div.fragment.well div.line:first').parent().removeClass('fragment well');
  }

  $('table.memberdecls').find('.memItemRight').each(function() {
    $(this).contents().appendTo($(this).siblings('.memItemLeft'));
    $(this).siblings('.memItemLeft').attr('align', 'left');
  });
	
  $(".memitem").removeClass('memitem');
  $(".memproto").removeClass('memproto');
  $(".memdoc").removeClass('memdoc');
  $("span.mlabel").removeClass('mlabel');
  $("table.memberdecls").removeClass('memberdecls');
  $("[class^=memitem]").removeClass('memitem');
  $("span.mlabels").removeClass('mlabels');
  $("table.mlabels").removeClass('mlabels');
  $("td.mlabels-right").removeClass('mlabels-right');
  $(".navpath").removeClass('navpath');
  $("li.navelem").removeClass('navelem');
  $("div.ah").removeClass('ah');
  $("div.header").removeClass("header");
  
  $('.mdescLeft').each(function() {
    if($(this).html()=="&nbsp;") {
      $(this).siblings('.mdescRight').attr('colspan', 2);
      $(this).remove();
    }
  });
  $('td.memItemLeft').each(function() {
    if($(this).siblings('.memItemRight').html()=="") {
      $(this).attr('colspan', 2);
      $(this).siblings('.memItemRight').remove();
    }
  });



  // replace directory icon
  $("span.iconfopen").addClass("glyphicon glyphicon-folder-open");
  $("span.iconfopen").removeClass("iconfopen");
  // replace file icon
  $("span.icondoc").addClass("glyphicon glyphicon-file");
  $("span.icondoc").removeClass("icondoc");

  // replace namespace and class label
  $("span.icona > span.icon").each(function() {
    $(this).addClass("label");
    $(this).removeClass("icon");
    $(this).after("&nbsp;");
    if($(this).text()=="N") $(this).addClass("label-warning");
    if($(this).text()=="C") $(this).addClass("label-info");
  });
  $("span.icona").removeClass("icona");

  // add icons to navbar1
  $('#navrow1 > ul > li > a[href="index.html"]').prepend('<span class="glyphicon glyphicon-home"></span>&nbsp;');
  $('#navrow1 > ul > li > a[href="pages.html"]').prepend('<span class="glyphicon glyphicon-paperclip"></span>&nbsp;');
  $('#navrow1 > ul > li > a[href="namespaces.html"]').prepend('<span class="label label-warning">N</span>&nbsp;');
  $('#navrow1 > ul > li > a[href="annotated.html"]').prepend('<span class="label label-info">C</span>&nbsp;');
  $('#navrow1 > ul > li > a[href="files.html"]').prepend('<span class="glyphicon glyphicon-duplicate"></span>&nbsp;');
  // add icons to navbar2 namespaces
  $('#navrow2 > ul > li > a[href="namespaces.html"]').prepend('<span class="glyphicon glyphicon-align-justify"></span>&nbsp;');
  $('#navrow2 > ul > li > a[href="namespacemembers.html"]').prepend('<span class="glyphicon glyphicon-list"></span>&nbsp;');
  // add icons to navbar2 classes
  $('#navrow2 > ul > li > a[href="annotated.html"]').prepend('<span class="glyphicon glyphicon-align-justify"></span>&nbsp;');
  $('#navrow2 > ul > li > a[href="classes.html"]').prepend('<span class="glyphicon glyphicon-sort-by-alphabet"></span>&nbsp;');
  $('#navrow2 > ul > li > a[href="hierarchy.html"]').prepend('<span class="glyphicon glyphicon-align-left"></span>&nbsp;');
  $('#navrow2 > ul > li > a[href="functions.html"]').prepend('<span class="glyphicon glyphicon-list"></span>&nbsp;');
  // add icons to navbar2 files
  $('#navrow2 > ul > li > a[href="files.html"]').prepend('<span class="glyphicon glyphicon-align-justify"></span>&nbsp;');
  // add icons to navbar3 namespaces/members
  $('#navrow3 > ul > li > a[href="namespacemembers.html"]').prepend('<span class="glyphicon glyphicon-th"></span>&nbsp;');
  $('#navrow3 > ul > li > a[href="namespacemembers_func.html"]').prepend('<span class="glyphicon glyphicon-option-horizontal"></span>&nbsp;');
  $('#navrow3 > ul > li > a[href="namespacemembers_type.html"]').prepend('<span class="glyphicon glyphicon-option-vertical"></span>&nbsp;');
  // add icons to navbar3 classes/members
  $('#navrow3 > ul > li > a[href="functions.html"]').prepend('<span class="glyphicon glyphicon-th"></span>&nbsp;');
  $('#navrow3 > ul > li > a[href="functions_func.html"]').prepend('<span class="glyphicon glyphicon-option-horizontal"></span>&nbsp;');
  $('#navrow3 > ul > li > a[href="functions_enum.html"]').prepend('<span class="glyphicon glyphicon-option-vertical"></span>&nbsp;');
});
