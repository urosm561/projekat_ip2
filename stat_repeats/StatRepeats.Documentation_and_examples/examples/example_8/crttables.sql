connect to ripiti;
drop index fragment__idx;
drop index sequence__idx;
drop table sequence;
drop table fragment;
drop table match;

create table sequence 
      (
       id          integer     not null generated always as identity ( start with 1 increment by 1 minvalue 1 maxvalue 2147483647 no cycle cache 20 no order ), 
       name        varchar(80) not null ,
       version     smallint    not null ,
       min_length  integer     not null ,  
       type        varchar(2)  not null , 
       finished    integer     not null , 
       confidence  decimal(4,3)                 ,
       primary     key(id)              ,
       check (type in ('dn','dc','in','ic'))

      ) not logged initially  compress yes;

create index sequence__idx on sequence ( name, version, min_length, type, confidence ) 
allow reverse scans page split symmetric collect sampled detailed statistics compress no;
      
comment on table sequence is 'The names and ordinal numbers of the genomes/proteins';      
comment on sequence (id         is 'Sequence id',
		     name       is 'Sequence name',
		     version    is 'Sequence version',
                     min_length is 'Minimal sequence length (parametar of the clled program)',
                     type       is 'Type od repeat: dc,dn,ic,in',
                     finished   is 'Reserve - direct insert into database',
                     confidence is 'Confidence of results obtained'
                    ); 

create table fragment 
      (
       id            integer not null generated always as identity ( start with 1 increment by 1 minvalue 1 maxvalue 2147483647 no cycle cache 20 no order ), 
       text          varchar(1000) not null, 
       text_length   integer                ,
       primary key   (id)
      ) not logged initially compress yes;

comment on table fragment is 'Sequences that are fragments of palindromes/repeats';      
comment on fragment (id          is 'Fragment id',
		     text        is 'Fragment sequence',
		     text_length is 'Fragment length'
                    );      
		     
create unique index fragment__idx on fragment ( text ) allow reverse scans page split symmetric collect sampled detailed statistics compress no;

       
create table match 
             (
              id             integer not null generated always as identity ( start with 1 increment by 1 minvalue 1 maxvalue 2147483647 no cycle cache 20 no order ), 
              seqid          integer not null , 
              id_fragment1   integer not null , 
              id_fragment2   integer not null , 
              location1      integer not null , 
              location2      integer not null ,
              primary key(id)                 ,
              foreign key(seqid)        references sequence,
              foreign key(id_fragment1) references fragment,              
              foreign key(id_fragment2) references fragment                            
             ) not logged initially compress yes;

comment on table match is 'Positions of found pairs';
comment on match (id           is 'Primary key of the relation',
		  seqid        is 'Id of processed genome/protein that contains fragment',
		  id_fragment1 is 'Id of left (first) component of repeat pair',
		  id_fragment2 is 'Id of right (second) component of repeat pair',
		  location1    is 'Position of left (pair) component in genome/protein',
		  location2    is 'Position of left (pair) component in genome/protein'
		  );      


RUNSTATS ON TABLE SEQUENCE ON ALL COLUMNS WITH DISTRIBUTION ON ALL COLUMNS AND SAMPLED DETAILED INDEXES ALL SET PROFILE;
RUNSTATS ON TABLE match ON ALL COLUMNS WITH DISTRIBUTION ON ALL COLUMNS AND SAMPLED DETAILED INDEXES ALL SET PROFILE;
RUNSTATS ON TABLE fragment ON ALL COLUMNS WITH DISTRIBUTION ON ALL COLUMNS AND SAMPLED DETAILED INDEXES ALL SET PROFILE;
commit;       
