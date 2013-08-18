/*global $,document,window*/
"use strict";

!function ($) {
    $("<style type='text/css'> .tree-node{ display:inline-block; width:1em; margin-right:0.5em; text-align: right;cursor: pointer;} " +
      ".tree-node-label{ display:inline-block; text-align: left;} </style>").appendTo("head");

    var Node = function (element, options) {
        this.element = $(element);

        // If the node has a subroot / children
        this.subroot = this.element.children('ul').root(options);

        // Value of the leaf node: if any
        this.text = this.element.clone().children('ul').remove().end().text();
        if (this.subroot === undefined)
            this.value = this.element.attr ('input-value') || this.text;


        // State of the current node
        this.state = $.fn.node.State.UNCHECKED;

        // Setup animation speed
        this.animation = this.element.data('tree-animation') || options.animation;


        // Type of node (checkbox? radio )
        this.type = options.type;


        if(this.type === 'checkbox') {
            // Create checkbox area
            var inputSpan = $('<span class="tree-node tree-checkbox" style=""></span>');

            // Create a checkbox and copy attributes from options
            var input = $('<input  class="node-input" type="checkbox"></input>');
            if(this.value)
                input.attr('value', this.value);

            $.each(options.input, function (key, value) {
                input.attr(key, value);

            })

                
                inputSpan.prepend(input);


            var custom = this.element.data('input-custom') || options.input.custom;
            if(custom != undefined && custom != undefined && custom != undefined && custom != undefined){
                this.custom = options.input.custom;

                var customSpan = $('<span class="node-custom"></span>');
                input.attr('type', 'hidden');

                inputSpan.prepend(customSpan);
            }

            this.element.prepend(inputSpan);


            if(this.custom){
                this.customElement = this.element.find('.node-custom').eq(0); 
                this.customElement.html(this.custom.unchecked);
            }

            this.checkbox = this.element.find('.tree-checkbox').eq(0);
            this.input = this.element.find('.node-input').eq(0); 

            this.checkbox.on({
                click: $.proxy(this.inputclick, this)
            });
        }


        if (this.subroot === undefined)
            this.element.prepend('<span class="tree-node tree-bullet"> <!--&#183;--> </span>');
        else 
            this.element.prepend('<span class="tree-node tree-bullet"><i class="icon-caret-down" ></i>  </span>');
        
        this.icon = this.element.find('.tree-bullet').eq(0);

        this.element.css('margin', '0');

        if (options.toggleLI == true){
            this.element.on({
                click: $.proxy(this.click, this)
            });
        }
        else {
            this.icon.on({
                click: $.proxy(this.click, this)
            });
        }


    }

    Node.prototype = {
        constructor: Node,

        isSubrootVisible: function () {
            var visible = false;
            if(this.subroot != null){
                this.subroot.each(function (){
                    if($(this).data('root').element.is(":visible"))
                        visible = true;
                });
            }
            return visible;
        },

        switch: function () {
            if(this.subroot != null){
                this.subroot.each(function () {
                    $(this).data('root').switch(this.animation);
                });
                this.icon.find('i').eq(0).toggleClass('icon-caret-right');
                this.icon.find('i').eq(0).toggleClass('icon-caret-down');
            }
        },
        collapse: function(recursive){
            if(this.subroot != null){
                this.subroot.each(function(){

                    $(this).data('root').collapse();
                    if(recursive == true){
                        $(this).data('root').collapseAll();
                    }
                });    

                this.icon.find('i').eq(0).addClass('icon-caret-right');
                this.icon.find('i').eq(0).removeClass('icon-caret-down');
            }
        },
        expand: function(recursive){
            if(this.subroot != null){
                this.subroot.each(function(){

                    $(this).data('root').expand();
                    if(recursive == true){
                        $(this).data('root').expandAll();
                    }
                });    
                this.icon.find('i').eq(0).addClass('icon-caret-right');
                this.icon.find('i').eq(0).removeClass('icon-caret-down');            
            }

        },

        update: function(state){
            if(state == $.fn.node.State.CHECKED){
                this.state = $.fn.node.State.CHECKED;
                this.input.prop('checked', true);
                this.input.prop('indeterminate', false);
                if(this.custom)
                    this.customElement.html(this.custom.checked);
            }
            if(state == $.fn.node.State.UNCHECKED){
                this.state = $.fn.node.State.UNCHECKED;
                this.input.removeAttr('checked');
                this.input.prop('indeterminate', false);
                if(this.custom)
                    this.customElement.html(this.custom.unchecked);
            }
            if(state == $.fn.node.State.INDETERMINED){
                this.state = $.fn.node.State.INDETERMINED;
                this.input.removeAttr('checked');
                this.input.prop('indeterminate', true);
                if(this.custom)
                    this.customElement.html(this.custom.indetermined);

            }

        },
        check: function(){

            this.update($.fn.node.State.CHECKED);
            // propagate downwards
            if(this.subroot){
                this.subroot.children('li').each(function(){
                    $(this).data('node').propagate($.fn.node.State.CHECKED);    
                });
            }
            // update parents
            var parent = this.element.parent().parent().data('node');
            if(parent){
                parent.propagate();
            }
        },
        uncheck: function(){
            
            this.update($.fn.node.State.UNCHECKED);

            // propagate downwards
            if(this.subroot){
                this.subroot.children('li').each(function(){
                    $(this).data('node').propagate($.fn.node.State.UNCHECKED);    
                });
            }
            // update parents
            var parent = this.element.parent().parent().data('node');
            if(parent){
                parent.propagate();
            }
            
        },

        inputclick: function(ev){
            //console.log(this);
            ev.stopPropagation();
            
            // unchecked and indetermined default to checked
            if(this.state == $.fn.node.State.CHECKED)
                this.uncheck();
            else 
                this.check();

            

        },

        propagate: function(state){
            if(state === undefined){

                var parent = this.element.parent().parent().data('node');
                if(parent)
                    parent.propagate();

                var inputs = this.element.find('input[value]');

                

                var a = $.map(inputs, function(val, i){
                    if($(val).prop('checked') === true)
                        return true;
                });

                if(a.length == inputs.length)
                    this.update($.fn.node.State.CHECKED);
                else if(a.length == 0)
                    this.update($.fn.node.State.UNCHECKED);
                else
                    this.update($.fn.node.State.INDETERMINED);

            }
            else{
                this.update(state);

                if(this.subroot){
                    this.subroot.children('li').each(function(){
                        $(this).data('node').propagate(state);    


                    });
                }

            }

        },

        click: function(ev){
            ev.stopPropagation();
            this.switch();    
        },

        search: function(str){
            var searchResult = false;

            if(str == "")
                searchResult = true;

            if(this.text.toLowerCase().indexOf(str.toLowerCase()) != -1)
                searchResult = true;


            if(this.subroot){
                this.subroot.children('li').each(function(){
                    if($(this).data('node').search(str))    
                        searchResult = true;
                });
            }

            if(searchResult){
                this.element.show(this.animation);
                if(this.subroot){
                    this.subroot.each(function(){
                        $(this).data('root').expandAll();
                    });
                }
            }
            else{
                this.element.hide(this.animation);
                if(this.subroot){
                    this.subroot.each(function(){
                        $(this).data('root').collapseAll();
                    });
                }

            }
            
            return searchResult;

        }

    };

    $.fn.node = function ( option, val ) {
        if(this.length == 0)
            return undefined;
        return this.each(function () {
            var $this = $(this),
            data = $this.data('node'),
            options = typeof option === 'object' && option;

            if (!data)  {
                $this.data('node', (data = new Node(this, $.extend({}, $.fn.node.defaults, options))));
            }
            if (typeof option == 'string') {
                data[option](val);
            }
        });
    };
    $.fn.node.defaults = {
    };

    $.fn.node.State = {
        CHECKED : 'checked',
        UNCHECKED : 'unchecked',
        INDETERMINED : 'indetermined'
    };

    $.fn.node.Constructor = Node;


    var Root = function(element, options) {
        this.element = $(element);

        this.element.css('list-style', 'none');
        this.element.css('margin-left', '0');
        this.element.css('padding-left', '1em');
        this.element.css('text-indent', '0em');
        this.element.css('display', 'block');
        
        this.nodes = this.element.children('li').node(options);

        this.animation = this.element.data('tree-animation') || options.animation;
    };

    Root.prototype = {
        constructor: Root,

        switch: function(){
            this.element.toggle(this.animation);
        },
        collapse: function(){
            this.element.hide(this.animation);
        },
        expand: function(){
            this.element.show(this.animation);
        },
        toggleAll: function(){
            var isVisible = false;
            if(this.nodes){
                $.each(this.nodes,function(){
                    if($(this).data('node').isSubrootVisible())
                        isVisible = true;
                });
            }
            if(isVisible)
                this.collapseAll();
            else
                this.expandAll();
        },
        collapseAll: function(){
            if(this.nodes){
                $.each(this.nodes,function(){
                    $(this).data('node').collapse(true);
                });
            }
        },
        expandAll: function(){
            if(this.nodes){
                $.each(this.nodes,function(){
                    $(this).data('node').expand(true);
                });
            }
        },
        selectAll: function(){
            if(this.nodes){
                $.each(this.nodes,function(){
                    $(this).data('node').check();
                });
            }
        },
        unselectAll: function(){
            if(this.nodes){
                $.each(this.nodes,function(){
                    $(this).data('node').uncheck();
                });
            }
        },
        search: function(str){
            //console.log(str);
            this.expandAll(0);
            if(this.nodes){
                var searchResult = false;
                $.each(this.nodes,function(){
                    if($(this).data('node').search(str))
                        searchResult = true;
                });
                return searchResult;
            }
        }
    }
    $.fn.root = function ( option, val ) {

        if(this.length == 0)
            return undefined ;
        return this.each(function () {
            var $this = $(this),
            data = $this.data('root'),
            options = typeof option === 'object' && option;
            if (!data)  {
                $this.data('root', (data = new Root(this, $.extend({}, $.fn.root.defaults, options))));
            }
            if (typeof option == 'string') {
                data[option](val);
            }
        });

    };

    $.fn.root.defaults = {
        type: 'none', // radio | checkbox | none 
        attrs: [],
        animation: 100,
        toggleLI: false,
    };

    $.fn.root.Constructor = Root;

    $.fn.tree = function ( option, val ) {

        if(this.length == 0)
            return undefined ;

        return $.each(this,function () {

            var $this = $(this);
            var data = $this.data('tree');
            var options = typeof option === 'object' && $.extend(true, {}, option);

            if(option === undefined) options = new Object;
            
            if (!data)  {


                if($this.is('input')){
                    $(this).wrap('<div></div>');
                    this.div = $(this).parent();
                }
                if($this.is('div')){
                    this.div = $this;
                }

                options.animation = $this.data('tree-animation') || options.animation || $.fn.root.defaults.animation;
                options.toggleLI = $this.data('tree-toggleli') || options.toggleLI || $.fn.root.defaults.toggleLI;
                options.collapsed = $this.data('tree-collapsed') || options.collapsed || $.fn.root.defaults.collapsed;
                options.selected = $this.data('tree-selected') || options.selected || $.fn.root.defaults.selected;

                options.type = $this.data('tree-type') || options.type || $.fn.root.defaults.type;
                if($this.is('input')){
                    options.type = 'checkbox';
                }


                options.attrs = $this.prop("attributes");

                // get input options from either data or options
                if(options.type === 'checkbox'){

                    if(options.input === undefined){

                        options.input = new Object();

                        $.each($this.prop("attributes"), function(){
                            if(this.name.substring(0, 15) == "data-tree-input"){
                                try{
                                    options.input[this.name.substring(16)] =  $.parseJSON(this.value);
                                }
                                catch(err){
                                    options.input[this.name.substring(16)] = this.value;
                                }
                            }
                        });

                        if($this.is('input')){
                            $.each($this.prop("attributes"), function(){
                                if(this.name.substring(0, 15) != "data-tree-input"){
                                    options.input[this.name] = this.value;
                                }
                            });
                        }
                    }

                }


                // Get the values for div and input there should be values to be 
                if($this.is('div') || $this.is('input')){
                    this.values = $this.data('tree-values') || options.values;
                    if(this.values === undefined){
                        $.getJSON( $this.data('tree-values-url')).done(function( data ) {/*console.log('should read json');*/} );
                    }

                    this.ul = $.fn.tree.generate(this.values, $.extend({},options));

                    this.div.append(this.ul);
                }

                // Wrap ul in div because why not
                if($this.is('ul')){
                    $(this).wrap('<div></div>');
                    this.ul = $this;
                }


                // CREATE the tree - roots and nodes and all...
                $this.data('tree', (data = this.ul.root($.extend({}, $.fn.root.defaults, options))));
                


                // Remove input tag which now is unnecessary
                if($this.is('input')){
                    $this.attr('type','hidden');
                }

                if(options.collapsed == true)
                    data.data('root').collapseAll(val);

                if(options.selected == true)
                    data.data('root').selectAll(val);
            }


            if (typeof option == 'string') {
                data.data('root')[option](val);
            }

        });

    }

    $.fn.tree.generate = function(data, options){

        if($(data['child']).length == 0)
            return;

        var $ul = $('<ul></ul>');

        $(data['child']).each(function(index, value){
            var $li = $('<li></li>');
            var $span = $('<span></span>');

            $span.html(value['text']);
            $span.addClass('tree-node-label');

            $li.attr('input-value', value['value']);
            $li.css('vertical-align','middle');

            $li.append($span);

            $li.append($.fn.tree.generate(this, $.extend({},options)));
            
            $ul.append($li);

        });

        return $ul;


    }

}( window.jQuery );
