dojo.provide("plugins.core.Template");

dojo.require("dojo.parser");


// REFERENCE SITES

// Object Oriented concepts and inheritance
// http://dojotoolkit.org/book/dojo-book-0-4/part-3-dojo-programming-model/object-oriented-concepts-and-inheritance

// BIOMART
//http://www.ensembl.org/biomart/martview/e61a39e4e8e306354f8f2a7a70bbc53c/e61a39e4e8e306354f8f2a7a70bbc53c

// DOJO DEMOS
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/demos/
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/dojox/form/tests/test_SelectStack.html

//Real World Dojo part 5: Custom Components
//http://www.dotnetmafia.com/blogs/jamesashley/archive/2008/10/28/761.aspx

dojo.declare(
    "plugins.core.Template",
    [dijit._Widget, dijit._Templated],
    {
        //Path to the template of this widget. dijit._Templated uses this to
        //snatch the template from the named file via a synchronous call.
        // This provides the relative path to the template file for the widget. Note that
        // fetching the template for a widget requires a synchronous network call, although
        // Dojo will cache the template string after it is fetched (unless a custom build
        // with interned template strings is used).
        templatePath: null,

        // Calls dijit._Templated.widgetsInTemplate
        // If widgets are defined inside of the template (path or string), this
        // value should be specified and set to true so that the Dojo parser will
        // find and instantiate those widgets. This value is false by default.
        // While this property is easy to overlook, know that defining widgets
        // inside of other widget?s templates can often be very useful.
        widgetsInTemplate : true,

        // SLOTS TO BE POPULATED BY THE DESCENDENT OF THIS CLASS
        title : '',
        id : '',
		filename: '',


		// STORE ELEMENTS OF INTEREST
		elements : [],
		elementHash: null,
        elementObjects : [],
		elementValues : null,

		// WIDGET TYPES FOR WHICH THE DOJO CLASSES ARE KNOWN
		knownTypes : null,

		// HASH OF FUNCTIONS TO GET AND SET VALUES IN WIDGETS AND DOM NODES
		valueFunctions : null,

		// CSS FILES
        cssFiles : [],

        //Any initialization code would go here in the constructor. dijit._Widget and
        //dijit._Templated do not have parameters in their constructors, so
        //there wouldn't be any multiple-inheritance complications
        //if you were to include some paramters here
        constructor : function (args)
        {

            this.id         = args.id;
			this.filename   = args.filename;
			//this.url        = args.url;

            // LOAD CSS FOR WIDGETS
            this.loadCSS();
        },


        //Inherited from dijit._Widget and called just before template
        //instantiation in buildRendering. This method is especially useful
        //for manipulating the template before it becomes visible.
        postMixInProperties: function()
        {
        },


        //You can override this method to manipulate widget once it is
        //placed in the UI, but be warned that any child widgets contained
        //in it may not be ready yet.        
        postCreate: function()
		{
		},


        //Called after the widget's children and all other widgets on the
        //page have been created. Provides an opportunity to manipulate child
        //widgets before they're displayed.
        // NOW THAT ALL THE CHILD WIDGETS HAVE BEEN CREATED, ADD EACH CHILD WIDGET
        // TO ITS RELATED widgetObject IN THE GLOBAL this.elementObjects HASH

		// USE this.inherited TO CALL THIS AT THE START OF THE SUBCLASS OBJECT'S
		// startup METHOD THEN SET UP THE PAGE WITH MORE COMPLICATED LOGIC 
        startup : function ()
        {

			// SET ELEMENT HASH
			this.setElementHash();

            // BUILD HASH OF ELEMENTS (name: ELEMENT OBJECT) AND SET
			// UNIQUE IDS FOR EACH 
            this.setElementObjects();

			// SET ELEMENTS: ADD SPECIFIC PROPERTIES, BASIC LOGIC AND
			// VALUES TO WIDGETS OR NODES
			this.setElements();
        },


		// setElements
		//
		// ADD SPECIFIC PROPERTIES AND VALUES TO WIDGETS OR NODES
		// 
		// USES THE FACT THAT this[functionName] CAN BE A FUNCTION

		setElements : function ()
		{

			for ( name in this.elementHash )
			{

				if ( this[name] )
				{
					// RUN THE FUNCTION
					this[name](this.elementObjects[name], name);

					// **** DEBUG : SHOW THE VALUE 
					if ( this.elementObjects[name] == null )
					{
						alert("this.elementObjects[" + name + "] == null");
					}

					if ( this.elementObjects[name].valueFunction == null )
					{
						this.elementObjects[name].valueFunction = this.valueFunction('null');
					}
				}
				else
				{
				}
			}
		},


		// SET HASH OF ELEMENTS OF INTEREST WITH IDS IN TEMPLATE.HTML
		setElementHash : function ()
		{

			// SET HASH OF ELEMENT IDS
			this.elementHash = new Object;
			for ( var i = 0; i < this.elements.length; i++ )
			{
				if ( this.elementHash[this.elements[i]] )
				{
				}
				else
				{
					this.elementHash[this.elements[i]] = 1;
				}
			}
		},


		// setElementObjects
		//
		// BUILD HASH OF ELEMENTS (name: ELEMENT OBJECT) AND SET
		//
		// UNIQUE IDS FOR EACH 
		//
		// NOTES:
		// dojo.query AND forEach
		// http://www.dojotoolkit.org/book/dojo-book-0-9/part-3-programmatic-dijit-and-dojo/functions-used-everywhere/dojo-foreach

        setElementObjects : function ()
        {
			var thisTemplate = this;

			dojo.query("*", this.containerNode).forEach(

				function(element) {

					if ( element.id && thisTemplate.elementHash[element.id] ) 
					{

						var widget = dijit.byId(element.id);
						if ( widget != null )
						{
							thisTemplate.elementObjects[element.id] = widget;
						}
						else
						{
							thisTemplate.elementObjects[element.id] = element;
						}


						// SET UNIQUE ID
						var uniqueId = dijit.getUniqueId(thisTemplate.id);
						element.id = uniqueId;
					}
				}
			);

			return;
        },




		// loadCSS
		//
		// LOAD CSSFILES IN this.cssFiles ARRAY

        loadCSS : function()
        {

            // LOAD CSS
            var cssFiles = this.cssFiles;

            for ( var i in cssFiles )
            {
                var cssFile = cssFiles[i];

                var cssNode = document.createElement('link');
                cssNode.type = 'text/css';
                cssNode.rel = 'stylesheet';
                cssNode.href = cssFile;
                cssNode.media = 'screen';
                document.getElementsByTagName("head")[0].appendChild(cssNode);
            }
        },




		// loadWidgetValues
		//
		// LOAD A LIST OF name : value HASHED VALUES INTO ELEMENTS IN THE TEMPLATE

        loadWidgetValues : function (response, ioArgs)
        {

			for ( var name in response )
			{

                var value = response[name];
				if ( this.elementObjects[name] )
				{
					this.elementObjects[name].value(this.elementObjects[name], value);
				}
				else
				{
				}
				//this.elementObjects[name].value(this.elementObjects[name], "'" + value + "'");
			}
        },

		// valueFunction
		//
		// GET OR SET THE VALUE OF A WIDGET OR NODE

		valueFunction : function (widget, name)
		{

			// SET VALUE FUNCTION FOR saveReport LATER
			switch (name)
			{
				case "combobox" : case "editor" : case "checkbox" :

					return function(value) {
						if ( value )
						{
							return widget.setValue(value);
						}
						else
						{
							var returnValue = widget.getValue();
							if ( returnValue != null && returnValue != '' )
							{
								returnValue = returnValue.replace(/"/g, "'");
							}
							return returnValue;
						}
					};

				case "null" :
					return function(value) { return null; };

				case "radio" : return function(value) {
						return value ? widget.checked = value : widget.checked;
					};

				case "textinput" : return function(value) {
						return value ? widget.value = value : widget.value;
					};

				case "div" : return function(value) {
						return value ? widget.innerHTML = value : widget.innerHTML;
					};			
			}

			return 0;
		}
    }


); // plugins.core.Template
