$(document).ready(function() {

    // Bootstrap tooltip settings
    $('body').tooltip({
        selector: '[data-toggle="tooltip"]',
        trigger: 'hover',
        placement: 'top',
        html: 'true'
    });

    // Bootstrap popover settings
    $('body').popover({
        selector: '[data-toggle="popover"]',
        trigger: 'hover',
        placement: 'right',
        html: 'true'
    });

    // bootstrap & django hack
    // alias for alert-error class -> alert-danger
    $('.alert-error').addClass('alert-danger');
})
