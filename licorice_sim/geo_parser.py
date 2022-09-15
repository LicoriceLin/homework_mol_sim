from typing import List,Dict
import torch
import re
from collections import OrderedDict

LENTHS={'angstrom':1.8897161646321,
        'bohr':1}

inner_empty_code=re.compile('^.*\s.+$')
def _check_unit(unit:str,unitdict):
    if unit in unitdict:
        pass
    else:
        raise ValueError(f'invalid unit. allows: {unitdict.keys()}')

def _to_unit(tensor,ori_unit,to_unit,units):
    scale=units[to_unit]/units[ori_unit]
    return tensor/scale

def split_(s:str):
    return s.split() if re.match(inner_empty_code,s) else [s.strip()]

class RawGeo:
    '''
    contains a OrderedDict,`self.segments`
    keys in the ordered dict should be List[List[str]]

    there must be some reserved keys for segments. To be update.
    '''
    def __init__(self) -> None:
        self.segments:Dict[str,List[List[str]]]=OrderedDict()

    @classmethod 
    def from_in(cls,in_file):
        cls=cls()
        key='holder'
        with open(in_file,'r') as f:
            for i in f.readlines():
                if i.startswith('%'):
                    key=i[1:-1]
                    assert key not in cls.segments.keys(),f'key redefined: {key}'
                    cls.segments[key]=[]
                else:
                    if len(i.strip())>0:
                        cls.segments[key].append(split_(i))
        return cls
    
    @classmethod
    def from_xyz(cls,xyz_file):
        cls=cls()
        with open(xyz_file,'r') as f:
            for i,j in enumerate(f.readlines()):
                if i==0:
                    cls.segments['ATOM_NUMBER']=[split_(j)]
                elif i==1:
                    cls.segments['META_INFO']=[split_(j)]
                elif i==2:
                    cls.segments['ATOMIC_POSTION']=[split_(j)]
                else:
                    cls.segments['ATOMIC_POSTION'].append(split_(j))
        return cls
    
    def __getitem__(self,__k)->List[List[str]]:
        return self.segments[__k]

    def out(self,outfile:str):
        if outfile.endswith('xyz'):
            self._out_xyz(outfile)
        elif outfile.endswith('in'):
            self._out_in(outfile)
        else:
            raise NotImplementedError

    def _out_xyz(self,outfile:str):
        '''
        a temporal method for hw1 
        to be rewrite.  
        '''
        with open(outfile,'w') as f:
            for key in ['ATOM_NUMBER','META_INFO','ATOMIC_POSTION','ATOMIC_POSTION']:
                if self.segments.get(key,None):
                    for l in self.segments[key]:
                        line='\t'.join(l)+'\n'
                        # import pdb
                        # pdb.set_trace()
                        f.write(line)

    def _out_in(self,outfile:str):
        with open(outfile,'w') as f:
            for k,v in self.segments.items():
                f.write(f'#{k}\n')
                for l in v:
                    line='\t'.join(l)+'\n'
                    f.write(line)


class Coordination:
    '''
    convert a raw text data of coord to tensor
    an extra para,unit is optional to handle unit-swtich issues. 
    (maybe a baseclass is better for all data-rich object?)
    '''
    def __init__(self,raw:List[List[str]],unit:str='angstrom') -> None:
        self.U=LENTHS
        _check_unit(unit,self.U)
        cord_str=[i[1:] for i in raw]
        cord_float = [ [ float(_) for _ in i] for i in cord_str]
        self.cord=torch.Tensor(cord_float)
        self.cord=self.cord.type(torch.float32)
        self.atom=[i[0] for i in raw]
        self.unit=unit

    def unit(self,unit:str)->None:
        _check_unit(unit,self.U)
        self.cord=_to_unit(self.cord,self.unit,unit,self.U)
        # scale=self.U[unit]/self.U[self.unit]
        # self.cord=self.cord/scale
        self.unit=unit

    def to_raw(self)->List[List[str]]:
        pass


class CellParameter:
    '''
    save the parameter for cell,
    unit is optional like `Coordination`
    '''
    def __init__(self,raw:List[List[str]],unit:str='angstrom') -> None:
        self.U=LENTHS
        _check_unit(unit,self.U)
        self.unit=unit
        self.base=OrderedDict()
        def raw2base(r):
            bl=[float(i) for i in r]
            bt=torch.Tensor(bl)
            return bt.type(torch.float32)  
        self.base['x']=raw2base(raw[0])
        self.base['y']=raw2base(raw[1])
        self.base['z']=raw2base(raw[2])
        self.fullbase=torch.cat([i.view(1,-1) for i in self.base.values()],dim=0)

    def unit(self,unit:str)->None:
        _check_unit(unit)

class hw1Universe:
    '''
    temporal class for system in hw1.
    '''
    def __init__(self,infile:str):
        if infile.endswith('.in'):
            self._raw=None

        pass
