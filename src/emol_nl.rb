class SimpleNameList
    def initialize( _filename = nil, _key_field = "$", _end_field = "$end", _septoken = ",", _comments_regexp = /#|!/ )
        @filename = _filename
        @keyfield = _key_field
        @endfield = _end_field
        @septoken = _septoken
        @comments = _comments_regexp
    end
    
    ######################################################
    def readline
        if ( !@filename ) then
            fp = STDIN
        else
            fp = open( @filename )
        end
        str = ""
        while line = fp.gets do
            line = line.chop
            s = line.split( @comments )[ 0 ]
            if ( !s.nil? ) then
                str += s
            end
        end
        return str
    end
    
    ######################################################
    def get_token( str, sep )
#        printf "get_token: str: \"%s\"\n", str
        if ( str.nil? ) then
            return nil, nil
        end
        pos_sep  = str.index( sep )
        if ( pos_sep ) then
            pos_psep = pos_sep + sep.length
#            printf "get_token: psep, ppsep: %d, %d\n", pos_sep, pos_psep
            str_car = str[ 0...pos_psep ]
            if ( pos_psep == str.length ) then
                str_cdr = nil
            else
                str_cdr = str[ pos_psep...str.length ]
            end
        else
            str_car = nil
            str_cdr = str
        end
#        printf "get_token: str_car, str_cdr: \"%s\", \"%s\"\n", str_car, str_cdr
        return str_car, str_cdr
    end

    ######################################################
    def query_fieldname( str )
        name = ( str.split )[ 0 ].downcase
        return name
    end

    ######################################################
    ######################################################
    def parse_field( str_field, item )
        name_field = query_fieldname( str_field )
        strcar = ""
        strcdr = str_field[ name_field.length...str_field.length ].strip
        name = ""
        while strcar do
#            printf "parse_field: cdr:      \"%s\"\n", strcdr
            strcar, strcdr = get_token( strcdr, @septoken )
            if ( strcar ) then
                strcar = strcar[ 0...(strcar.length - @septoken.length ) ].strip
                if ( strcar.index( "=" ) ) then
                    data = []
                    strs = strcar.split( "=" )
                    name = strs[ 0 ].downcase
                    if ( !item[ name ] ) then
                        item[ name ] =  []
                        item[ name ] << strs[ 1 ].strip
                    else
                        item[ name ] << strs[ 0 ].strip
                    end
                else
                    if ( name.length > 0 ) && ( item[ name ] ) then
                        item[ name ] << strcar
                    end
                end
            end
            if ( strcdr ) then
                strcdr = strcdr.strip
#                printf "parse_field: car, cdr: \"%s\", \"%s\"\n", strcar, strcdr
            end
        end
    end

    ######################################################
    def parse_allfield( fields )
        items = Hash.new
        fields.each do |field|
            name = query_fieldname( field )
            if ( !items[ name ] ) then
                items[ name ] = Hash.new
            end
            parse_field( field, items[ name ] )
        end
        return items
    end

    ######################################################
    def split_to_field( str )
        fields = []
        pre_substr = str[ 0...str.length ].split( @keyfield )
        for i in 1...pre_substr.length
            fields << pre_substr[ i ]
        end
        return fields
    end

    ######################################################
    def parse
        strfile = readline 
        fields = split_to_field( strfile )
        items  = parse_allfield( fields )
        return items
    end

end
