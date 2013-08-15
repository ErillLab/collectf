
var block_message = '<h1>Please wait.<br/><i class="icon-spinner icon-spin"></i></h1>';



$(document).ready(function() {
    $(document).on('click', 'a.block_before_load', function(ev) {
	ev.preventDefault();
	console.log($(this).attr('href'));
	$.blockUI({message: block_message});
	window.location = $(this).attr('href');
    });

});
